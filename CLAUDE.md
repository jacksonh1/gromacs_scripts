# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Project overview

Two GROMACS pipelines — **T-REMD** and **plain production MD** — for characterizing
**designed protein structures** on single-node GPU clusters (SLURM). Developed for the
Keating lab at MIT; configurable for any cluster via `site_config.sh`.

**Every run starts from a folded/designed input structure — never from an unfolded or
extended state.** The tools serve four objectives:

1. Stability/rigidity characterization of a given starting (designed) structure
2. Identifying flexible regions (per-residue)
3. Variant comparison — run the same protocol on several design variants and see which best retains its designed conformation
4. Bound-state sampling — simulate a complex in its bound pose and sample the bound ensemble

**Tool mapping:** plain MD is mainly for **#4** (and optionally #1, #2); T-REMD is
primarily for **#1–3**. Neither is a folding-from-unfolded tool — the input pose is the
reference the analysis is measured against (RMSD = drift from the design, RMSF = local
flexibility).

**Entry point:** copy `example/submit_jobs/submit_REMD.sh` (T-REMD) or
`example/submit_jobs/submit_MD.sh` (plain MD), set parameters, run it. The submit script
submits the matching engine (`scripts/simulation/REMD-gromacs.sbatch` /
`MD-gromacs.sbatch`) to SLURM with parameters exported as environment variables.

**Configuration:** `site_config.sh` in the repo root — cluster-specific settings (GROMACS path, scratch root, module names). Edit once per cluster; all scripts source it automatically.

**Post-analysis:** automatically run at the end of each job (STEP 12 of the sbatch). Scripts live in `scripts/analysis/`; see `scripts/analysis/README.md` for usage. Reference for GROMACS T-REMD log format and parsing: `scripts/analysis/REMD_log_reference.md`.

**Example job:** `example/outputs/output_T-REMD/helix_fusion-2ns-REMD-300-400K-48reps-NVT-exf-1ps/` — use this to test analysis scripts without rerunning a simulation.

### Critical T-REMD concept

Each `prod/rep{i}/` directory is a **fixed temperature slot**, not a physical configuration. Coordinates exchange between slots at each exchange step. This means:

- `rep000/remd.xtc` is the 300 K constant-temperature ensemble — use it directly for analysis.
- **No demux step is needed** to obtain the thermodynamic ensemble at a given temperature.
- `gmx demux` follows a *configuration* as it walks through temperature space — useful for visualizing the random walk, but not for thermodynamic analysis.

### Pipeline stage folders are named by role, not ensemble

Output stages are named for what they *do*, so the name stays correct regardless of ensemble:

- **MD:** `build/ → em/ → heat/ → density/ → [relax/] → prod/`
- **REMD:** `build/ → em/ → density/ → equil/ → prod/`

(`heat/` = NVT thermalization, MD-only; `density/` = iterative NPT density equilibration; `relax/` = optional unrestrained NPT, MD-only, `RELAX_NS>0`; `equil/` = per-replica equilibration, REMD-only.) Do **not** reintroduce ensemble-named folders (`nvt/`, `npt/`, `eq/`). Production basenames stay `md.*`/`remd.*` (method, not ensemble) — the analysis layer keys off them.

### REMD ensemble toggle (`ENSEMBLE=NVT|NPT`)

T-REMD production runs NVT (default, constant volume) or NPT (`ENSEMBLE=NPT`, C-rescale barostat at `REF_P` bar). Under NPT the per-replica `equil/` stage *and* production are pressure-coupled, and GROMACS adds the *PV* term to the replica-exchange Metropolis criterion automatically — `-replex` is unchanged. The constant-temperature ensemble interpretation above is unaffected. Plain MD production is always NPT.

---

## Working style

### Push back when something doesn't make sense

If a request is scientifically or technically wrong — the wrong GROMACS flag, a misunderstanding of T-REMD, a workflow that would silently produce bad results — say so clearly before implementing it. A wrong simulation or analysis run wastes real compute time and can produce results that look plausible but are wrong. See the **Known gotchas** section below for specific discovered pitfalls.

### Suggest better alternatives

If there's a cleaner, faster, or more correct approach than what was asked for, say so and explain the tradeoff. Don't just implement what was asked if something clearly better exists. Include enough context for an informed decision.

### Prefer simple and scalable solutions

Solutions should work for 8, 48, or 128 replicas without special-casing. Prefer shell/Python idioms that stay readable as the codebase grows. If a task has a five-line solution and a fifty-line solution, understand why the complexity is or isn't justified before recommending it. Don't add abstractions or generality beyond what the current task requires.

---

## Code style: fail loudly

This project follows a "fail loudly" philosophy. Bugs that crash immediately are strongly preferred over bugs that silently produce wrong results.

### Core philosophy

- **Crashes are cheap; silent bugs are expensive.** Prefer code that crashes obviously when assumptions are violated over code that "handles" the violation by producing degraded output.
- **Don't paper over uncertainty.** If you're unsure whether a value can be None, empty, or wrong-typed, either ask, add an assertion, or leave a clearly-marked comment — never add a default to make the question go away.
- **Make illegal states unrepresentable.** Prefer types and structures where the invalid case can't be expressed, over runtime checks for the invalid case.

### Error handling

- **No bare `except:` or `except Exception:` clauses** unless the exception is logged AND re-raised, or the recovery is documented and intentional.
- **Don't catch exceptions just to log and continue.** If the operation failed, the caller needs to know.
- **No `.get(key, default)` patterns** unless the default is semantically meaningful, not just a way to avoid a `KeyError`.
- **No `value or fallback` shortcuts** (`x or []`, `x or {}`, `x or 0`) unless `None`/empty/zero is genuinely interchangeable with the fallback. These hide bugs where `x` was unexpectedly empty.
- **Don't add defensive `if x is not None:` checks** unless `None` is a real expected case. If `None` would be a bug, let it crash.

### Indexing and iteration

- **Prefer iteration over indexing.** Use `for item in items`, not `for i in range(len(items))`. When you need the index too, use `enumerate`.
- **Use `zip(strict=True)`** (Python 3.10+) so mismatched-length iterables crash instead of silently truncating.
- **Assert invariants before code that relies on them.** E.g., `assert len(a) == len(b)` before zipping when the lengths must match.

### Types and structure

- **Use dataclasses or TypedDicts, not raw dicts**, when the shape matters and is fixed.
- **Parse, don't validate, at boundaries.** Convert untrusted input into typed structures at the edge; the rest of the code should be able to assume the data is valid.

### When in doubt

- **Ask before adding error handling.** If tempted to wrap something in try/except, ask what the intended behavior is when it fails.
- **Flag assumptions explicitly.** If making an assumption about input shape, range, or type that isn't enforced by the types, leave a comment like `# ASSUMES: items is non-empty`.

---

## Known gotchas

Specific pitfalls discovered in this environment. **When a new pitfall is discovered during work — a surprising GROMACS behavior, a cluster quirk, a wrong assumption that caused a failure — add it here immediately without waiting to be asked.** The goal is that the same mistake is never made twice.

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

`gmx rms`, `gmx gyrate`, and `gmx rmsf` all write distances in **nm**; structural work expects ångström. `plot_xvg.py` converts any axis whose `.xvg` label carries an `nm` unit to Å (×10, label rewritten), leaving the underlying `.xvg` data in nm. If you add a new distance metric, confirm its plot reads in Å.

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

### GROMACS: conformational clustering uses sklearn (`cluster_traj.py`), not `gmx cluster`

`gmx cluster` (gromos) builds the full pairwise-RMSD matrix → **O(N²)** in time and memory, which is impractical for the long production runs (25k+ frames). Conformational clustering is therefore done in `scripts/analysis/cluster_traj.py` (MDAnalysis + scikit-learn DBSCAN/k-means on flattened Cα coordinates), which scales: a 25k-frame run clusters in ~13 s. It runs in the **shared** part of `run_analysis.sh` (consumes `<prefix>_stripped_aligned.{xtc,gro}`), so it serves single- and multi-chain alike — no `multichain_*` variant.

Things to keep straight:
- **`--cutoff` is a real RMSD cutoff (nm)** only because the input frames are pre-aligned to one common reference, so flattened-coord Euclidean distance = `√N_atoms × RMSD`. The script sets DBSCAN `eps = cutoff_Å × √N_selected`. **If you ever feed `cluster_traj.py` an un-aligned trajectory, the cutoff stops meaning RMSD.** (This was a latent bug in the Amber `cluster_MD.py`, where `eps` was a raw flattened distance mislabelled "Å".)
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

**Fix:** keep **`REPLEX_PS ≥ 1.0`** (1–2 ps is standard and sub-ps buys no extra sampling). Only if sub-ps exchange were truly required, add **`-update cpu`** to the STEP 9 `mdrun` (moves the swap off the GPU) at ~10–25% throughput cost. `config_example.sh` ships `0.5` — raise it.

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

