# Known gotchas — detailed reference

Pitfalls discovered in this environment, with the full explanation and fix for each. The
compact **index** (one line per gotcha, phrased as the trigger) lives in `CLAUDE.md` under
"Known gotchas" — that's the always-loaded tripwire; this file is the detail you open when a
tripwire fires. **When you discover a new pitfall: append the detail here AND add a one-line
entry to the index in `CLAUDE.md`** — both, every time, or the tripwire is invisible / the
detail is lost. The goal is that the same mistake is never made twice.

---

### SLURM: sbatch scripts are copied to a temp path at execution time

SLURM copies the `.sbatch` script to a temporary location before running it, so `BASH_SOURCE[0]` inside an sbatch script does not point to the repo. Code like `$(dirname "${BASH_SOURCE[0]}")/../site_config.sh` will silently resolve to the wrong place or fail.

**Fix (two options):**
- Use `$SLURM_SUBMIT_DIR` — a SLURM-provided environment variable that always holds the directory from which `sbatch` was called. Works as long as the job is submitted from the repo root (the normal case).
- Or pass the repo path explicitly via `--export` (e.g. `GROMACS_SCRIPTS_DIR`) and use that variable inside the sbatch.

Regular shell scripts called from *outside* SLURM (e.g. `scripts/analysis/`) can use `BASH_SOURCE[0]` reliably.

### GROMACS: `-pbc nojump` breaks T-REMD trajectories

`nojump` works by comparing consecutive frames to detect box-boundary crossings. T-REMD coordinate exchanges cause discontinuous jumps between frames that `nojump` misinterprets and incorrectly "fixes."

**Fix:** use `-pbc mol` instead. It is a per-frame operation and is safe for REMD trajectories.

### GROMACS: `gmx trjconv -s` must match the XTC atom count exactly

If the structure passed to `-s` has fewer atoms than the XTC (e.g. a protein-only GRO when the XTC has the full solvated system), GROMACS silently truncates the read, produces a degenerate fitting matrix, and dies with `Too many iterations in routine JACOBI`.

**Fix:** always use a full-system GRO (matching the XTC) as `-s`. Extract a separate protein-only GRO if needed for downstream analysis.

### GROMACS: a CUDA "illegal memory access" (#700) is usually a physics blow-up, not a GPU fault

When a fully GPU-resident run (`-nb gpu -pme gpu -bonded gpu`, single rank → GPU-resident update) dies at step 0 with `CUDA error #700 (cudaErrorIllegalAddress)` — often surfacing as an assertion in `freeDeviceBuffer`/`cudaStreamQuery` — the GPU is rarely the real problem. A numerical explosion (overlapping atoms → NaN coordinates → out-of-bounds pairlist) manifests as an illegal memory access on the GPU instead of a clean error.

**Fix:** re-run the *same tpr* with `-nb cpu`. The CPU path prints the real diagnosis — e.g. `Constraint error in algorithm Lincs at step 0` plus an energy table showing a huge **positive** `LJ (SR)` and absurd temperature/pressure (clashes). Then fix the physics (usually a bad starting structure), not the GPU flags. In one case the culprit was the heat-stage `grompp` reading the *un-minimized* `build/*_ions.gro` instead of `em/em.gro`, discarding EM and starting on clashing coordinates. Always confirm each equilibration `grompp -c` consumes the *previous* stage's output (MD: em.gro → heat.gro → density_seg*.gro → [relax.gro] → prod; REMD: em.gro → density_seg*.gro → equil.gro → prod), not the build.

### GROMACS: a fixed `SEED` reproduces the setup, NOT a bit-for-bit GPU trajectory

The pipeline takes a `SEED` job parameter (default `-1` = random per run, i.e. legacy behaviour; `>= 0` pins it). When set it fills `gen-seed` (velocity generation), `ld-seed` (the **stochastic** `V-rescale` thermostat *and* `C-rescale` barostat draw from it — not just Langevin/SD integrators), and `mdrun -reseed` (the replica-exchange RNG, counter-based on `(seed, step)`). Per-replica stages offset the base as `SEED + i` (via the `rep_seed` helper) so replicas stay decorrelated; single-sim stages (MD `heat`, `density`) use `SEED` directly. `REMD-restart.sbatch` recovers `-reseed` from the original `mdrun_remd.log` so a seeded run resumes reproducibly.

**Do not promise bitwise reproducibility.** These are GPU-resident runs (`-nb gpu -pme gpu -bonded gpu`); GROMACS does floating-point reductions in a nondeterministic order on GPU (and under domain decomposition), so two same-seed runs on the same hardware still diverge within a few ps — MD is chaotic. What a fixed seed *does* guarantee: identical initial velocities and identical RNG streams, i.e. reproducible **initial conditions and setup** (useful for variant comparison starting from the same thermal kick) and full provenance in `parameters.txt`. Bitwise determinism would require `-nb cpu -pme cpu` single-rank — far too slow to be worth it here. Adding a seed to a new mdp block means BOTH `gen-seed` (if `gen-vel = yes`) and `ld-seed` (any V-rescale/C-rescale block); `SEED` is a shared param, so it lives in all three `validate_params.py` lists + each engine's STEP-1 read/export (see the STEP-1 validation gotcha).

### Bash: `bash -n` does not catch `set -u` unbound variables — grep removed vars when refactoring

The sbatch engines run `set -euo pipefail`. A reference to a variable that was never assigned (e.g. left dangling after you delete the block that defined it) is a **runtime** crash — `line N: VAR: unbound variable` — that `bash -n` (syntax-only) passes clean. This bit us once: refactoring the MD post-analysis step removed `ANALYSIS_DIR=...`, but the STEP 11 summary still echoed `${ANALYSIS_DIR}`, so the job ran to completion and then died on the last summary line.

**Fix:** when you remove or rename a variable's assignment, `grep -n 'VARNAME' the_script` and confirm there are no remaining references (or none beyond the assignment). Don't rely on `bash -n`; it won't flag this. (Note the failure surfaces *after* the real work — the analysis had already run — so a late summary crash does not mean the outputs are missing.)

### GROMACS: `-pbc mol`/`-pbc whole` split multi-chain complexes — multi-chain has its own pipeline

For a multi-chain complex, `-pbc mol` wraps each chain's COM into the box *independently*, and `-pbc whole` only un-breaks within a molecule — neither keeps the two chains in the *same* periodic image. So a bound complex can end up with one chain a box-vector away ("split"), which: inflates the whole-complex RMSD against a split reference (one case: median backbone RMSD 24.5 Å, while Rg stayed ~1.4 nm — the giveaway that it's a relative-image artifact, not real motion), and biases Rg/RMSF. Whether it splits is frame-dependent, so the trajectory flickers between intact and split.

**Architecture:** the **single-chain** analysis path is left simple — `fix_PBC.sh` uses `-pbc mol` and `extract_protein.sh` uses `-pbc whole`. **Multi-chain** systems are handled by a **separate `multichain_*` pipeline** that `run_analysis.sh` selects automatically when the topology has >1 protein chain. Those scripts apply the **`-pbc cluster`** periodic-image fix (cluster group = `Protein`), which pulls all chains into one image and is per-frame → REMD-safe (unlike `nojump`). **Note:** `-pbc cluster` is a PBC operation; it is unrelated to *conformational* clustering (`gmx cluster`) — the scripts are named `multichain_*`, never `*cluster*`, to avoid that confusion.

**Caveat:** `-pbc cluster` assumes the complex stays within ~half the (minimum) box vector. Beyond that the periodic image is genuinely ambiguous *and* the box is too small (minimum-image violation → corrupt physics), so it is a setup red flag, not just an analysis issue. The multi-chain pipeline emits an inter-chain minimum-distance curve (`gmx mindist`, minimum-image-aware) so you can confirm the complex stayed bound; for genuine dissociation studies use a bigger box.

### GROMACS: `gmx rms` does not make its reference whole — a PBC-broken reference inflates RMSD

`gmx rms -s REF -f TRAJ` least-squares-fits *every frame* onto REF before measuring, so the reference does **not** need to be pre-aligned to the trajectory (the fit is internal and per-frame). But `gmx rms` does **not** apply any PBC treatment to REF. If the reference structure (e.g. a protein extracted straight from `em.gro`) has the molecule split across a box boundary, the reference coordinates are physically wrong and *every* RMSD value is silently corrupted — typically showing an absurd, near-constant offset. In one example this inflated the t=0 backbone RMSD to ~13 Å; with a whole reference the same trajectory read a physical ~3 Å.

**Fix:** make the reference whole when extracting it — `trjconv -s system.tpr -pbc whole`. `extract_protein.sh` does this. It is a no-op if the molecule was already intact, so it is safe to always apply. (Note: a near-constant, physically-impossible RMSD that does *not* start at ~0 is the signature of this bug, not of real drift.)

### GROMACS: distances are nm — convert to Å for analysis plots

`gmx rms`, `gmx gyrate`, and `gmx rmsf` all write distances in **nm**; structural work expects ångström. `gromd-plot-xvg` converts any axis whose `.xvg` label carries an `nm` unit to Å (×10, label rewritten), leaving the underlying `.xvg` data in nm. If you add a new distance metric, confirm its plot reads in Å.

### GROMACS: `gmx demux` is not for extracting constant-T ensembles

`gmx demux` reconstructs the trajectory of a single *configuration* as it walks through temperature space. It does not produce the thermodynamic ensemble at a given temperature.

**Fix:** use the slot trajectory directly. `prod/rep000/remd.xtc` is already the 300 K constant-temperature ensemble — no post-processing needed to obtain it.

### GROMACS: REMD acceptance rates are pre-computed in the log

The `Replica exchange statistics` block at the end of each replica log contains pre-computed per-pair acceptance rates, exchange counts, and mean Metropolis probabilities. Reparsing the thousands of per-frame `Repl ex` lines to recount exchanges is unnecessary.

**Fix:** parse the statistics block at the end of the log. See `scripts/analysis/REMD_log_reference.md` for the format and parsing code.

### GROMACS: the Empirical Transition Matrix is not dwell time

The matrix printed after the statistics block records one-step transition *probabilities* between temperature slots. It is not a histogram of time spent at each temperature.

### GROMACS: this build provides only `gmx_mpi` — there is no `gmx` or `gmx mdrun`

The install (`$HOME/opt/gromacs/2024.3-plumed/bin/`) ships a single binary, `gmx_mpi`. There is no plain `gmx`, so `gmx mdrun` does not exist. All commands — preprocessing *and* mdrun — must use `gmx_mpi`. An MPI-compiled GROMACS runs serial steps fine as a single rank with no `mpirun` (this is how the REMD engine runs EM/density); only multi-replica steps need `mpirun -np N`.

**Fix:** use the `GMX="${GMX:-gmx_mpi}"` pattern and call `$GMX mdrun` everywhere. (Both engines previously carried a dead `MDRUN="gmx mdrun"` variable that was never used and would fail here — it has been removed; don't reintroduce it.) The analysis scripts already probe `gmx_mpi` first, then fall back to `gmx`.

### GROMACS: conformational clustering uses sklearn (`gromd-cluster`), not `gmx cluster`

`gmx cluster` (gromos) builds the full pairwise-RMSD matrix → **O(N²)** in time and memory, which is impractical for the long production runs (25k+ frames). Conformational clustering is therefore done in `gromd_analysis/clustering.py` (the `gromd-cluster` entry point) (MDAnalysis + scikit-learn DBSCAN/k-means on flattened Cα coordinates), which scales: a 25k-frame run clusters in ~13 s. It runs in the **shared** part of `run_analysis.sh` (consumes `<prefix>_stripped_aligned.{xtc,gro}`), so it serves single- and multi-chain alike — no `multichain_*` variant.

Things to keep straight:
- **`--cutoff` is a real RMSD cutoff (nm)** only because the input frames are pre-aligned to one common reference, so flattened-coord Euclidean distance = `√N_atoms × RMSD`. The script sets DBSCAN `eps = cutoff_Å × √N_selected`. **If you ever feed `gromd-cluster` an un-aligned trajectory, the cutoff stops meaning RMSD.** (This was a latent bug in the Amber `cluster_MD.py`, where `eps` was a raw flattened distance mislabelled "Å".)
- **Cluster count vs noise is controlled by `min_samples`** (the DBSCAN density knob), **not** a post-hoc population filter — the user rejected adding one. The default is **adaptive: `max(10, 1.5% of frames)`**, so a region must hold ~1.5% of the trajectory to be a state and the long tail of tiny clusters falls into noise (without it, a flexible system gave 80 clusters). Raise `--min-samples` for fewer; pass an absolute int to override. Tuned on the WW-domain REMD slots to keep even the most heterogeneous case to ≲10 clusters.
- **Outputs go in a `clustering/` subdir** next to the prefix (`analysis/clustering/<prefix>_cluster_*`), not flat in `analysis/`.
- **`-pbc cluster`** (the multi-chain periodic-image fix) is **unrelated** to this conformational clustering — different operation, despite the shared word. The multichain PBC scripts are deliberately not named `*cluster*`.

### Per-job parameters are validated at STEP 1 (fail loud, not silent default)

Both engines call `scripts/simulation/validate_params.py` (`--engine remd|md`) in STEP 1: once with `--check-keys "$1"` *before* sourcing a config file (rejects typo'd keys against an allowlist), then once on the resolved values (required `PDB_IN`, type, range — covers both config-file and `--export` params). A bad parameter now exits with `[ERROR] …` instead of silently using a `${VAR:-default}`. `PDB_IN` is **required** (no default).

Things to keep straight when touching parameters:
- **Three lists must stay in sync** when you add/rename a param: the `${VAR:-default}` read in the engine's STEP 1, the `export …` list right before the `python3 "$VALIDATE_PARAMS"` call (so the validator can read it from the environment), and the allowlist + rules in `validate_params.py`. A param missing from the export list is silently unvalidated; one missing from the allowlist is falsely rejected as a typo.
- **`validate_params.py` is stdlib-only** — it runs under the sbatch's *system* `python3` at STEP 1, before the analysis conda env is activated. Don't add non-stdlib imports.
- **Derived step counts** (`*_NSTEPS`, `TAU_T`, `REPLEX_STEPS`) are computed **after** the validator on purpose, so a non-numeric input fails with a clear message instead of a raw traceback from the inline `python3 -c`. Keep new derivations below the validate call.
- **Typo'd-key detection is config-file only.** Under `--export=ALL` the job inherits the whole shell environment, so there's no clean manifest to flag an unknown `--export` key against — only value/required/range checks cover that path (by design).

### GROMACS: short `REPLEX_PS` (<1 ps) deadlocks GPU-resident REMD — keep it ≥ 1 ps

On GPU-resident REMD production (`-nb/-pme/-bonded gpu`, GPU update, many ranks/GPU — e.g. 48 ranks on 4 L40S), `REPLEX_PS=0.5` (exchange every 250 steps @ 2 fs) reliably **hangs mid-run**: CPUs at 100%, GPUs at 0%, indefinitely. Signature: *all* replicas frozen at the *same* exchange step with **clean physics** (`Constr. rmsd = 0`, no LINCS/NaN) — an **MPI collective deadlock at exchange**, not a blow-up (contrast the CUDA #700 gotcha). Cause: each accepted swap moves GPU-resident coordinates via GPU↔CPU transfers + stream syncs; at 250-step spacing consecutive swaps' stream work overlaps and races. It's a **threshold, not linear**: the *same* config at `REPLEX_PS=1.0` never hangs — interval was the only change.

**Fix:** keep **`REPLEX_PS ≥ 1.0`** (1–2 ps is standard and sub-ps buys no extra sampling). Only if sub-ps exchange were truly required, add **`-update cpu`** to the STEP 9 `mdrun` (moves the swap off the GPU) at ~10–25% throughput cost. **Both engine defaults are `1.0`** (REMD's was raised from `0.5` on 2026-09-02, along with `scripts/simulation/config_example.sh`), so a job now has to lower it deliberately.

### SLURM/GROMACS: `gmx_mpi` needs `libhcoll`/`libocoms` — missing on some `mit_normal_gpu` nodes; also constrain to `rocky8`

The `gmx_mpi` binary (`~/opt/gromacs/2024.3-plumed`) carries a **baked-in `NEEDED` on `libhcoll.so.1` and `libocoms.so.0`** (from `/opt/mellanox/hcoll/lib`), left over from its build-time MPI. These are **inert at runtime** — the loaded `openmpi/5.0.8` module is built `--without-hcoll`, so MPI actually routes through `libmpi.so.40` — but the loader still has to *resolve* them at exec. `/opt/mellanox` is **node-local**, and some `mit_normal_gpu` nodes (seen on node4207) don't have it, so the job dies **instantly** with `gmx_mpi: error while loading shared libraries: libhcoll.so.1: cannot open shared object file`. On `pi_keating` (and the login node) it resolves via `/etc/ld.so.conf.d/hcoll.conf`, which is why REMD there never hit this.

**Fix (two parts):**
- **Stage the two libs on shared storage and prepend `LD_LIBRARY_PATH`** in the job, *after* sourcing GMXRC: copy `libhcoll.so.1{,.0.9}` + `libocoms.so.0{,.0.0}` from `/opt/mellanox/hcoll/lib` to a shared dir, then `export LD_LIBRARY_PATH="<that dir>:${LD_LIBRARY_PATH:-}"`. Their transitive deps are all standard RDMA/system libs (`libibverbs`, `librdmacm`, `libnl`, …) present on any IB node, so no rabbit hole.
- **Add `#SBATCH -C rocky8`.** The binary is built on Rocky 8 (login node = Rocky 8.10, el8); landing on a non-rocky8 node would be a real glibc/ABI mismatch. All `mit_normal_gpu` nodes are rocky8 today, so this is future-proofing, not a live bug. (Discovered while smoke-testing the REST2 `-hrex` engine on `mit_normal_gpu`; see `knowledgebase/plans/REST2-pipeline.md`.)

### PLUMED/GROMACS: `-hrex` (REST2 Hamiltonian exchange) is SILENTLY BROKEN on the `2024.3-plumed_2.9.4` build

The current build (**GROMACS 2024.3 + PLUMED 2.9.4**) accepts `-hrex`, patches cleanly, loads PLUMED ("GROMACS-like replica exchange is on"), and *runs to completion* — but **every exchange has `dE_term = 0.000e+00` and zero replicas ever swap.** So a REST2 run would look successful while actually being N independent MD runs with no enhanced sampling — a silent-wrong result. Verified 2026-07-29 on both CPU and GPU with correctly-scaled `partial_tempering` topologies (the scaling itself is fine; `scale=1.0` reproduces the FF to the bit). Root cause: the PLUMED **2.9.x** hrex patch for **GROMACS 2024.x** is an incomplete port — `md.cpp`'s `GREX cacheLocalUSwap` energy re-evaluation returns a wrong value on the 2024 force API, so `replicaexchange.cpp`'s `if(plumed_hrex) delta=0.0` is never corrected by a real Hamiltonian delta. Confirmed authoritatively by PLUMED developer **Giovanni Bussi** on the plumed-users list: **"hrex is not yet supported in the 2024 patch"** — i.e. the `gromacs-2024.3` patch exists and supports general PLUMED use (biasing/metadynamics), but the hrex energy machinery was never ported to it. Corroborated by GROMACS forum + PLUMED GitHub issue #1326 (open as of Oct 2025, no fix). Per Bussi, **GROMACS 2023 is the newest GROMACS with working PLUMED hrex** — this is not an arbitrary downgrade. The only confirmed-working combo is **GROMACS 2023.5 + PLUMED 2.9.x** (scale-1.0 sanity check gives P=1.0 there); GROMACS 2025 + PLUMED 2.10 has a new *native* interface but is unverified for hrex.

**Fix:** do NOT ship a REST2 pipeline against this build. Rebuild a GROMACS+PLUMED with a working hrex patch first, and **gate the REST2 engine on an acceptance-rate self-check** (assert nonzero exchanges in a short pre-run, fail loud) so a regression can never silently produce a fake REST2. Plain T-REMD and MD on `2024.3-plumed_2.9.4` are unaffected (they don't use `-hrex`). See `knowledgebase/plans/REST2-pipeline.md`.

### PLUMED: the `plumed` CLI needs the hcoll libs too (same fix as `gmx_mpi`)

The standalone `plumed` CLI (`plumed info`, `plumed patch`, `plumed partial_tempering`) links **the exact same `libhcoll.so.1` / `libocoms.so.0`** as `gmx_mpi`, so it dies the same way on GPU nodes lacking `/opt/mellanox/hcoll`: `plumed: error while loading shared libraries: libhcoll.so.1`. (Confusingly, the plumed *kernel* `libplumedKernel.so` loads fine at `mdrun -plumed` time on those nodes — it's only the CLI binary that's affected.)

**Fix:** the same `LD_LIBRARY_PATH="$HCOLL_COMPAT_DIR:…"` staging that fixes `gmx_mpi`. The REST2 engine sets it before STEP 6, so `plumed partial_tempering` runs inline on the compute node — no need to pre-generate topologies elsewhere. **Exception — the GROMACS build:** `install_gromacs-2023.5-plumed.sh` does NOT set the hcoll path during the compile job, so it instead requires the source to be **pre-patched on the login node** (`plumed patch -p -e gromacs-<ver>`); the runtime/"shared" patch is only needed at patch time, not compile time, so cmake/make/install then run fine on the compute node against the pre-patched source (the script checks for `.preplumed` backups and errors with instructions otherwise).

### GROMACS/SLURM: compile on the same (or older) CPU microarch as the run nodes — SIMD/`-march` is baked into the binary

GROMACS auto-detects the **build node's** CPU at configure time and bakes its SIMD level + compiler `-march` into the binary (our builds: `-march=skylake-avx512`, `SIMD = AVX_512`). Run that binary on a node whose CPU **lacks** an instruction it was compiled with and it dies with **"illegal instruction"** (or GROMACS refuses at startup on a SIMD mismatch). This cluster is heterogeneous: the `pi_keating` build nodes **node[3619-3620]** are Skylake-AVX512, while `mit_normal_gpu` mixes L40S/H100/H200 hosts, some on newer microarchs (e.g. Sapphire Rapids). Building on a newer node and running on `node[3619-3620]` is the dangerous direction.

**Fix:** build on **node[3619-3620]** (`-p pi_keating --nodelist=node[3619-3620]`), OR build anywhere and force the target with `-DGMX_SIMD=AVX_512 -DCMAKE_C_FLAGS=-march=skylake-avx512 -DCMAKE_CXX_FLAGS=-march=skylake-avx512` (what `install_gromacs-2023.5-plumed.sh` does, so it can build on a free node yet stay portable).

### GROMACS: `OUTDIR` must be on the same filesystem as the submit directory (solvate/genion topology rename)

`gmx solvate` / `gmx genion` update the topology by writing a temp file (`temp.topXXXXXX`) **in the current working directory** and then `rename()`-ing it to the `-p` path. `rename(2)` cannot cross filesystems, so if `OUTDIR/build/` is on a *different* mount than the job's CWD (= `SLURM_SUBMIT_DIR`), the build dies at STEP 3 with a cryptic `System I/O error: Failed to rename temp.top... to .../<name>.top` (and leaves the topology missing, only a `#...top.1#` backup). This bit the REST2 engine test when `OUTDIR` was set to `/orcd/data/keating/...` while submitting from `/orcd/pool/004/...` — two different mounts. Affects all engines (REMD/MD/REST2), not just REST2.

**Fix:** keep `OUTDIR` on the **same filesystem** as where you `sbatch` from (the normal case — `OUTDIR` under the repo, submitted from the repo). Trajectories still go to `SCRATCH_ROOT` on any filesystem (those are written directly / via symlink, not renamed). A Skylake-AVX512 binary runs on both those nodes *and* the newer `mit_normal_gpu` nodes (a newer CPU is a superset), so building on the *oldest* target microarch is portable everywhere; the reverse is not. If you must build on a newer node but run on Skylake, pin `-DGMX_SIMD=AVX_512` **and** the compiler arch (`-DCMAKE_C_FLAGS=-march=skylake-avx512`, same for CXX). (Surfaced building the REST2 GROMACS 2023.5; the 2024.3 build already lived on node[3619-3620].)

### Bash/GROMACS: `VAR=$(cmd | grep -c ...)` under `set -e` silently kills the job when `cmd` fails to launch

The REST2 engine's STEP 0 hrex probe was `_HREX_OK=$("$GMX" mdrun -h 2>&1 | grep -ci 'hrex')`. Two failure modes compound: (1) if `gmx_mpi` **can't run at all** on the compute node (bad MPI/CUDA/library env), its real error is redirected by `2>&1` **into the pipe** and swallowed by `grep`; (2) `grep -c` then finds zero matches and **exits 1**, which — because it's the command substitution in an assignment — trips `set -e` and exits the whole script. Net result: the job dies in ~6 s having printed **only** `===== JOB START =====`, with an **empty** stderr and exit code 1 — a maximally uninformative failure (the `set +o pipefail` around it guards SIGPIPE from `grep -q`, but does nothing for the zero-match case). Seen 2026-07-29: a REST2 test on `pi_keating` node1927 failed exactly this way; the same binary/env ran clean on other rocky8 nodes (16-task allocation, GPUs) — i.e. the underlying `gmx_mpi` fault was **node-specific/transient**, but the engine hid it completely.

**Fix (applied):** never let a launch failure hide inside a `| grep` under `set -e`. Run the probe once, capturing output **and** exit code with `set -e` disabled around it, check the exit code explicitly, and **print the binary's real output** before exiting if it's nonzero — only then grep the captured text for the feature. Pattern:
```bash
set +e; _HELP="$("$GMX" mdrun -h 2>&1)"; _RC=$?; set -e
(( _RC == 0 )) || { echo "[ERROR] '$GMX mdrun -h' exited $_RC — binary could not run:"; echo "$_HELP"; exit 1; }
grep -qi 'hrex' <<<"$_HELP" || { echo "[ERROR] no -hrex in this build"; exit 1; }
```
General rule: a `VAR=$(... | grep -c/-q ...)` whose upstream command can fail to execute is a silent-death trap under `set -e`; separate "did it run?" from "what did it say?".

### Output model: folder-symlink stages — don't symlink `build/`, and don't copy the symlinks

The engines write the **bulk stage dirs** (`density/ equil/ prod/`, plus MD `heat/ [relax/]`) as **folder symlinks** `${OUTDIR}/<stage> → ${SCRATCH_DIR}/<stage>`, pre-created *before* mdrun so trajectories land on scratch without touching the pool (default `SYMLINK_BULK=1`; `SYMLINK_BULK=0` makes them real in `OUTDIR`). The transparency is the whole point — every stage-dir variable still reads `${OUTDIR}/<stage>`, so mdrun/grompp/analysis paths are unchanged. Three traps to keep in mind:

- **Never symlink `build/`** (or `em/`). `build/` must be **real in `OUTDIR`** or the `solvate`/`genion` topology `rename()` becomes cross-filesystem and dies (see the OUTDIR gotcha above). It's small; keep it real.
- **The stage symlink must PRE-EXIST before mdrun**, created in the tree-setup block. If a later step does `mkdir -p "$STAGE_DIR"` on a path that isn't already the symlink, it creates a *real* dir in OUTDIR and the bulk lands on the pool. (MD's optional `relax/` is created the same conditional way inside STEP 6.)
- **`finalize_outputs.sh` must skip symlinks.** It copies the real OUTDIR dirs into `SCRATCH_DIR` to make the archive self-contained; it explicitly `continue`s on any entry that is a symlink (`[[ -L ]]`), because copying `prod/` (which points back into `SCRATCH_DIR`) would recurse/duplicate. A plain `cp -a prod scratch/` would also nest as `scratch/prod/prod`.

**Fix / rule:** when touching the output tree, preserve these invariants — real `build/ em/ analysis/ logs/` (+ REST2 `topol/`), symlinked bulk stages created up front, and a symlink-skipping finalize. Recovering `SCRATCH_DIR` from a finished run (e.g. the restart engine): `SCRATCH_DIR="$(dirname "$(readlink -f "${OUTDIR}/prod")")"` when `prod/` is a symlink.

### GROMACS: CHARMM force fields REQUIRE force-switched vdW — plain-cutoff mdp silently gives wrong energetics

The engines' mdp nonbonded block is AMBER-style: `rvdw = rcoulomb = CUTOFF_NM` (1.0 nm) plain Verlet cutoff + `DispCorr = EnerPres`. CHARMM36/36m was **parameterized with force-switched van der Waals** (`vdw-modifier = force-switch`, `rvdw-switch = 1.0`, `rvdw = rcoulomb = 1.2 nm`) and **no analytic dispersion correction** (`DispCorr = no` — the force-switch already handles the vdW tail; adding DispCorr double-counts). Running CHARMM with the AMBER cutoff settings *runs fine and looks plausible* but produces subtly wrong forces/energies — the classic silent-bad-result failure. Equally, pairing CHARMM protein with plain TIP3P (or AMBER protein with the CHARMM-modified TIP3P) is wrong: each force field is validated only with its matched water. In the GROMACS CHARMM port, `pdb2gmx -water tip3p` resolves to the **CHARMM-modified TIP3P automatically** (it reads `<ff>.ff/tip3p.itp`, i.e. the LJ-on-H variant), so water follows the force field for free — never force plain TIP3P onto CHARMM.

**Fix (applied):** all three engines gate the nonbonded mdp lines on `[[ $FF == charmm* ]]` — CHARMM emits the force-switch block + `DispCorr=no` and forces `CUTOFF_NM=1.2` (with an `[INFO]` override notice); the AMBER branch (`VDW_BLOCK=""`) leaves the generated mdp **byte-identical** to before, so existing amber runs are unchanged. The CHARMM36m force field itself is not bundled with GROMACS (only `charmm27.ff` is) — it's the MacKerell **force-switch** port (`charmm36-feb2026_cgenff-5.0`) installed under `$HOME/opt/gromacs/ff/` and found via `export GMXLIB` in `site_config.sh` (searched in addition to each build's bundled `top/`, so amber keeps working). **Do NOT use the `ljpme` port** with these force-switch settings — LJ-PME needs `vdwtype = PME` + a different mdp path. Select per job with `FF=charmm36-feb2026_cgenff-5.0` (keep `WATER=tip3p`). Verified: pdb2gmx builds the topology ("The Charmm36m force field and the tip3p water model are used") and grompp accepts the force-switch tpr.

**Keep the upstream dated directory name.** The port was briefly installed here as `charmm36m.ff`, which is *not* how these ports are distributed — every release ships a dated directory (`charmm36-jul2022.ff`, `charmm36-feb2026_cgenff-5.0.ff`, …) and you pass that name to `pdb2gmx -ff`. Renaming to a generic `charmm36m` reads better in a submit script but makes `parameters.txt` record `FF=charmm36m` forever, so two runs a year apart log an identical string while a newer port has been dropped into the same directory — an invisible force-field change in a pipeline whose purpose is comparing variants. The mdp gate is a `charmm*` glob, so the dated name needs no code change. (The port's own provenance survives inside `forcefield.itp`'s header — charmm2gmx version + build date — but that describes the directory's *current* contents, not what a past job used.) For usability, `site_config.sh` maps short names to installed ports (`FF_ALIASES`, e.g. `charmm36m`); the engines resolve the alias at STEP 1 via `resolve_ff` and record the **resolved dated name** in `parameters.txt`, so the short name is input sugar only and never reaches the run record. An alias whose target is not in `GMXLIB` fails at STEP 1.

### GROMACS: `gmx trjconv` cannot combine `-fit` with `-pbc mol` — they must be separate passes

`trjconv` refuses `-fit` together with any `-pbc` mode other than `none`/`whole`, because a PBC wrap applied after a rotation would re-wrap atoms against a box the coordinates no longer match. This is *why* the analysis layer is split into `fix_PBC.sh` (PBC) and `strip_and_align_trajectory.sh` (fit) rather than one script — the two-script shape is a GROMACS constraint, not a stylistic choice.

**Fix / rule:** any new tool that needs both PBC treatment and alignment must run two `trjconv` calls, PBC first. See `extract_solvated_snapshots.sh` for the pattern.

### GROMACS: `-fit rot+trans` rotates coordinates but NOT the box — never run a PBC-aware selection on a fitted frame

This is the subtle one. `gmx trjconv -fit rot+trans` rewrites the coordinates but leaves the box vectors **byte-identical**. The frame is therefore internally inconsistent: PBC-aware distance code still measures minimum images against the *original* box, which no longer corresponds to the rotated coordinates. Any `gmx select` distance keyword (`within`, `same residue as (within …)`) run on a fitted frame silently returns a **wrong** atom set.

Measured on `example/outputs/output_MD/helix_fusion-2ns-MD-300K-NPT` at t=500 ps with a 0.5 nm shell: selecting on the fitted frame gave **1622** atoms vs **1598** on the same frame pre-fit — 8 spurious waters, the furthest **23 Å** from the protein, which render as water floating in empty space. Nothing errors; the PDB just quietly contains the wrong molecules.

**Fix / rule:** compute the selection on the **PBC-corrected but unfitted** frame, then apply the resulting index to the fitted frame — `-fit` does not renumber atoms, so the index stays valid. Order is **PBC → select → fit → write**. More generally: treat any fitted/rotated frame as having a meaningless box, and never feed one to a distance-based selection or to `gmx mindist`/`gmx select`.

### GROMACS: a solvent-shell selection must come AFTER PBC centering

`gmx select`'s `within` is PBC-aware, so on a raw frame it will select a water whose *periodic image* is near the protein while its stored coordinate sits across the box — written out as a water far from everything. Running `-pbc mol -center -ur compact` first puts the protein at the box centre with all nearest images already in place, so the selected waters are the ones physically drawn next to it.

**Fix / rule:** PBC-correct (and for a complex, `whole → cluster → mol+center+compact`) before any shell selection. Sanity-check the result by measuring the max water→protein distance: it should not exceed the shell by more than ~1 Å (whole-residue selection means a water's O can sit slightly beyond a cutoff satisfied by its H).

### GROMACS: `gmx select` names a selection with a LEADING QUOTED STRING, not `Name =`

`-select 'Shell = group "Protein" or ...'` does **not** name a selection — in the GROMACS selection language `name = expr` declares a reusable *variable*, which is not itself a selection, so the command dies with `Error in user input: Too few selections provided`. The naming syntax is a leading quoted string: `-select '"Shell" group "Protein" or ...'`.

Second trap in the same area: for a **dynamic** selection, `gmx select -on` stamps the frame and time into the group name — `[ Shell ]` is written as `[ Shell_f0_t1000.000 ]`. Selecting that group by name downstream fails. Select it by **index** instead (`printf "0\n" | gmx trjconv ... -n shell.ndx`); one selection over one frame means the ndx holds exactly one group, so `0` is unambiguous.

### GROMACS: do NOT parse `gmx check` for a trajectory's last frame time

Two independent traps, which is why the answer is "use something else entirely":

1. **The output uses carriage returns.** `gmx check -f traj.xtc` draws progress with `\r`, so every `Reading frame … time …` update *and* any final `Last frame  200 time 2000.000` share **one physical line**. `awk '/^Last frame/'` never matches and silently yields an empty string.
2. **`Last frame` is throttled and often never printed at all.** The progress interval widens as the trajectory gets longer (every frame → every 10 → …). On a long run the final update is simply skipped. Measured on a 2501-frame, 808 MB xtc: the last thing printed was `Reading frame 2000 time 40000.000`, no `Last frame` line was emitted, and **`gmx check` still exited 0**. The same parse worked on 14 other jobs in the same sweep — i.e. it fails for *some* inputs only, the worst failure mode there is.

**Fix / rule:** dump the final frame and read the time GROMACS stamps into the `.gro` **title line**, which is written unconditionally:

```bash
printf "System\n" | $GMX trjconv -s "$TPR" -f "$XTC" -o last.gro -dump 999999999
LAST_PS=$(head -1 last.gro | sed -n 's/.*t=[[:space:]]*\([0-9.eE+-]\{1,\}\).*/\1/p')
```
(`gmx check`'s `Item / #frames / Timestep` table is a second option, but `(n-1)×dt` assumes a uniform timestep and a `t=0` start, so it is strictly weaker.) Costs one pass over the trajectory — the same as `gmx check` would.

More generally: **never parse a GROMACS progress readout.** It is throttled, carriage-return drawn, and goes to stderr. Parse a file GROMACS wrote.

### Bash: a `[[ … ]] && echo` as the LAST line of a `set -e` script fails the job

Every engine sbatch runs under `set -euo pipefail` and ends with a block of summary
`echo`s. Writing a conditional one of those as `[[ -d "$DIR" ]] && echo "..."` is a trap:
when the test is false the compound command returns 1, and because it is the *last*
command the script exits 1 — SLURM then reports a fully successful run as **FAILED**.
(`set -e` does not fire on the test itself, since it is the left side of `&&`; the damage
is entirely in the exit status.)

**Fix / rule:** use a real `if` for conditional output near the end of a script, or append
`|| true`. The same applies to any `((counter++))` (returns 1 when the pre-increment value
is 0) or `grep -q` used as the final statement.

### Python packaging: setuptools flat-layout auto-discovery fails in this repo

`gromd_analysis/` sits at the repo root (flat layout, not `src/`). setuptools' automatic
package discovery refuses to guess when the root holds several candidate directories, and
this repo has `scripts/ example/ docs/ knowledgebase/ logs/` alongside the package. Without
an explicit declaration `pip install -e .` dies at build time with:

```
error: Multiple top-level packages discovered in a flat-layout: [...]
```

**Fix:** `pyproject.toml` declares the package explicitly — do not remove this stanza, and
add to it (rather than deleting it) if a second package is ever added:

```toml
[tool.setuptools]
packages = ["gromd_analysis"]
```

Flat layout was chosen deliberately over `src/`: it makes `PYTHONPATH=<repo root>` a working
no-install fallback on a node without network. `src/` would need `PYTHONPATH=<repo>/src`,
which is easier to get wrong. The tradeoff `src/` normally buys — guaranteeing tests run
against the *installed* package rather than the working copy — is worth nothing here, since
the install is always editable and never a published wheel.

Related: install with `pip install --no-deps -e .`. conda already provides
matplotlib/numpy/MDAnalysis/scikit-learn from `environment.yml`; letting pip resolve them
shadows the conda builds with wheels.

### GROMACS: there is no `pdb2gmx -ff list` — it reads `list` as a force-field name

`gmx_mpi pdb2gmx -f x.pdb -ff list` does **not** print the available force fields (that
form exists in some other tools' docs and in older tutorials). GROMACS 2024 looks for a
force field literally named `list` and dies:

```
Could not find force field 'list' in current directory, install tree or GMXLIB path.
```

— followed by an `MPI_ABORT` on an MPI build, which looks alarming for what is a typo-class
error. **Fix:** run `pdb2gmx` with **no** `-ff` and read the interactive menu, or list the
directories directly: `ls "$GMXLIB"/*.ff` plus `<build>/share/gromacs/top/*.ff`. The real
check that a force field resolves is a throwaway `pdb2gmx` run with `-ff <name>`; on
success it prints `Using the <Name> force field in directory <path>`.

### PLUMED/REST2: `partial_tempering` does NOT scale CHARMM CMAP — CHARMM + REST2 is silently wrong

`plumed partial_tempering` scales the solute Hamiltonian by rewriting the processed topology.
Measured on a real CHARMM36m topology (helix_fusion, 433 marked solute atoms), diffing the
λ=1.0 against the λ=0.5 topology, the sections it changes are:

```
151434  [ pairtypes ]
   865  [ dihedrals ]
   565  [ atomtypes ]
   433  [ atoms ]
   333  [ nonbond_params ]
```

`[ cmap ]` (41 entries) and `[ cmaptypes ]` (1475 entries) are **byte-identical between λ=1.0
and λ=0.5** — the script has no `cmap` handling at all (`grep -i cmap partial_tempering.sh`
returns nothing). CMAP is CHARMM's backbone φ/ψ cross-term correction and applies **only to the
protein**, i.e. exactly the hot region, so a CHARMM REST2 run scales charges, LJ and dihedrals
on the solute while leaving its backbone conformational term at full strength. The solute
Hamiltonian is only partially scaled, so the effective-temperature ladder is not the one the
run reports — and the term left unscaled is the one governing secondary-structure sampling,
which is usually the whole point of the run.

**Neither existing self-check catches it.** The λ=1.0 replica is still exact, so the
`scale=1.0 → P=1.0` sanity pair passes; exchanges still happen at nonzero rates, so the hrex
acceptance gate passes too. It looks like a healthy REST2 run.

AMBER force fields in this pipeline (`amber99sb-ildn`, `amber14sb`) have no CMAP, so the
existing REST2 path is unaffected. (Note `ff19SB` *does* use CMAP — same trap if it is ever
added.)

**Fix:** treat `FF=charmm*` as unsupported for REST2 and fail at STEP 1 rather than produce a
plausible-looking wrong ensemble. T-REMD and plain MD have no such restriction — they scale
nothing, so CHARMM is fully correct there. Lifting the restriction means teaching
`partial_tempering` to scale `[ cmap ]`/`[ cmaptypes ]` (or pre-scaling the CMAP grids per
replica), then re-validating.

Unrelated to the GROMACS build: the 2023.5 REST2 build finds `GMXLIB` and builds a CHARMM36m
topology fine (verified). The problem is the topology *scaling* step, not force-field access.

**Generalize this before adding any force field.** CHARMM is one instance of a class:
`partial_tempering` scales `[ atoms ] [ atomtypes ] [ nonbond_params ] [ pairtypes ] [ dihedrals ]`
and silently passes through everything else, so *any* force field with an extra conformational
term is affected — CMAP (CHARMM, AMBER ff19SB), polarization (Drude), tabulated torsions.
`docs/FORCE_FIELDS.md` carries the rule and a GPU-free per-force-field check (diff the λ=1.0 and
λ=0.5 topologies by section, then justify every unchanged section). Run it once per new force
field; the `charmm*` guard catches only the case we knew about.
