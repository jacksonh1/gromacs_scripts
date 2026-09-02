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

**Re-running analysis on a finished job:** `bash scripts/analysis/run_analysis.sh OUTDIR [REP]` — one command for MD/T-REMD/REST2 (auto-detected from the layout), idempotent. For a long trajectory, submit it instead: copy `example/submit_jobs/submit_analysis.sh` (→ `scripts/analysis/analysis.sbatch`, CPU-only, no GPU requested). Knobs (`CLUSTER_CUTOFF`, `SHELL_NM`, `N_SNAPSHOTS`) are **env overrides, not job parameters** — deliberately, so `validate_params.py` and its three lists stay out of it.

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

### Output model: folder-symlinks into scratch (`SYMLINK_BULK`)

Default (`SYMLINK_BULK=1`): the **small, laptop-worthy** dirs (`build/ em/ analysis/ logs/`, `solvated_snapshots/` (MD only), `parameters.txt`, final PDB) are **real** in `OUTDIR`; the **bulk stage dirs** (`density/ equil/ prod/` for REMD/REST2; `heat/ density/ [relax/] prod/` for MD; `topol/` stays real for REST2) are **folder symlinks into `SCRATCH_DIR`**, created *before* the run so mdrun writes trajectories/intermediates straight onto scratch. This keeps bulk off the tight-quota pool and makes a laptop rsync of `OUTDIR` drag only a few broken folder-links instead of thousands of per-file links. **Key invariants:** every stage-dir variable (`PROD_DIR`, `DENSITY_DIR`, …) still points at `${OUTDIR}/<stage>` — the symlink is transparent, so **no `deffnm`/`-x`/`grompp` path changes** are needed and analysis reads `prod/…/…​.xtc` through the symlink. `build/` must stay **real** (the `solvate`/`genion` topology `rename()` must stay same-filesystem — see the OUTDIR gotcha). On success a shared `finalize_outputs.sh` copies the real OUTDIR dirs into `SCRATCH_DIR` so the scratch archive is self-contained. There is **no `trajectories/` collection dir** — it was dropped. `SYMLINK_BULK=0` makes every stage dir **real in `OUTDIR`** with no scratch offload (then `SCRATCH_DIR` is empty and all scratch-aware steps skip). When adding a param, keep the three lists in sync (see the STEP-1 validation gotcha); `SYMLINK_BULK` is a shared param (allowlist + `flag01`).

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

## Tests

`pytest` from the repo root (in `groMD_env`). `tests/` covers the **Python** half only —
the parsers in `gromd_analysis/`: `.xvg`, topology sections, `JobDir.detect` (all three
engines + failure modes), the replica-exchange statistics block, DSSP. Fixtures are
synthesized in `tmp_path`; **no test may read `example/outputs/`** — it is `.gitignore`d and
its bulk stage dirs are scratch symlinks that get purged (a shipped example directory
disappeared mid-session once). When you touch a parser, add or update its test.

The `.sh` steps are not unit-tested — verify them by running `run_analysis.sh` on a finished
job and diffing the outputs against the previous run. That is the regression check for
anything touching the GROMACS path.

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

Specific pitfalls discovered in this environment. **The full explanation + fix for each lives in
`knowledgebase/GOTCHAS.md`** — the list below is a one-line index (phrased as the trigger) so the
tripwire stays visible in every session; open the doc when one applies. **When a new pitfall is
discovered during work — a surprising GROMACS behavior, a cluster quirk, a wrong assumption that
caused a failure — add its detail to `GOTCHAS.md` immediately AND add a one-line entry here**
(both, every time). The goal is that the same mistake is never made twice.

- **SLURM: sbatch scripts are copied to a temp path** — `BASH_SOURCE[0]` inside an sbatch does *not* point to the repo; use `$SLURM_SUBMIT_DIR` or pass the path via `--export`.
- **GROMACS: `-pbc nojump` breaks T-REMD trajectories** — exchanges look like jumps; use `-pbc mol` (per-frame, REMD-safe).
- **GROMACS: `gmx trjconv -s` must match the XTC atom count** — a too-small `-s` GRO silently truncates → `JACOBI` crash; always pass a full-system GRO.
- **GROMACS: a CUDA "illegal memory access" (#700) is usually a physics blow-up, not a GPU fault** — re-run the same tpr with `-nb cpu` to see the real error (usually a clashing starting structure / wrong `grompp -c` input).
- **GROMACS: a fixed `SEED` reproduces the setup, not a bit-for-bit GPU trajectory** — `SEED` (default -1 = random) pins `gen-seed`/`ld-seed` (V-rescale *and* C-rescale are stochastic) + `mdrun -reseed`, per-replica offset `SEED+i`; GPU float-reduction order still makes trajectories diverge, so never promise bitwise determinism.
- **Bash: `bash -n` does not catch `set -u` unbound variables** — after deleting a var's assignment, `grep` the script for remaining references; a dangling `${VAR}` crashes at runtime (often in a late summary line, *after* the real work).
- **GROMACS: `-pbc mol`/`-pbc whole` split multi-chain complexes** — neither keeps chains in the same image; multi-chain has a separate `multichain_*` pipeline using `-pbc cluster` (auto-selected). `-pbc cluster` ≠ conformational clustering.
- **GROMACS: `gmx rms` does not make its reference whole** — a PBC-broken reference silently inflates *every* RMSD (near-constant offset, doesn't start at ~0); extract the ref with `-pbc whole`.
- **GROMACS: distances are nm — convert to Å for analysis plots** — `gmx rms`/`gyrate`/`rmsf` write nm; `gromd-plot-xvg` converts nm-labelled axes to Å. Confirm any new distance metric reads in Å.
- **GROMACS: `gmx demux` is not for extracting constant-T ensembles** — it follows a configuration through temperature space; use the slot trajectory (`rep000/remd.xtc`) directly.
- **GROMACS: REMD acceptance rates are pre-computed in the log** — parse the `Replica exchange statistics` block at the end; don't recount per-frame `Repl ex` lines.
- **GROMACS: the Empirical Transition Matrix is not dwell time** — it's one-step transition probabilities between slots, not time-per-temperature.
- **GROMACS: this build provides only `gmx_mpi`** — there is no `gmx`/`gmx mdrun`; use the `GMX="${GMX:-gmx_mpi}"` pattern and `$GMX mdrun` everywhere.
- **GROMACS: conformational clustering uses sklearn (`gromd-cluster`), not `gmx cluster`** (O(N²)) — `--cutoff` means RMSD only because inputs are pre-aligned; `min_samples` (adaptive) controls cluster-vs-noise.
- **Per-job parameters are validated at STEP 1 (`validate_params.py`)** — three lists must stay in sync when adding/renaming a param (engine default read, `export` list, allowlist+rules); stdlib-only; derive step counts *after* the validator.
- **GROMACS: short `REPLEX_PS` (<1 ps) deadlocks GPU-resident REMD** — MPI collective hang at exchange with clean physics; keep `REPLEX_PS ≥ 1.0` (or add `-update cpu`); all engine defaults are now 1.0, so a job must lower it deliberately.
- **SLURM/GROMACS: `gmx_mpi` needs `libhcoll`/`libocoms`** — missing on some `mit_normal_gpu` nodes → instant "cannot open shared object" death; stage the libs + prepend `LD_LIBRARY_PATH`, and add `#SBATCH -C rocky8`.
- **PLUMED/GROMACS: `-hrex` (REST2) is SILENTLY BROKEN on the `2024.3-plumed_2.9.4` build** — patches/runs but every `dE_term = 0`, zero swaps; needs GROMACS 2023.5 + PLUMED. Gate any REST2 run on an acceptance self-check.
- **PLUMED: the `plumed` CLI needs the hcoll libs too** — same `LD_LIBRARY_PATH` fix as `gmx_mpi` (the kernel loads fine at mdrun time; only the CLI binary is affected).
- **GROMACS/SLURM: compile on the same (or older) CPU microarch as the run nodes** — SIMD/`-march` is baked in; build on `node[3619-3620]` (Skylake-AVX512) or pin `-DGMX_SIMD=AVX_512` + `-march=skylake-avx512`.
- **GROMACS: `OUTDIR` must be on the same filesystem as the submit directory** — `solvate`/`genion` do a cross-fs-fatal topology `rename()`; keep `OUTDIR` under the repo you `sbatch` from.
- **Bash/GROMACS: `VAR=$(cmd | grep -c ...)` under `set -e` silently kills the job when `cmd` fails to launch** — separate "did it run?" (capture rc with `set +e`) from "what did it say?" (then grep).
- **Output model: folder-symlink stages** — bulk stage dirs (`prod/ equil/ density/` …) are symlinks into scratch (default `SYMLINK_BULK=1`); never symlink `build/`/`em/` (cross-fs rename), pre-create the stage symlinks before mdrun, and `finalize_outputs.sh` must skip symlinks (`[[ -L ]]`) or it recurses. See the "Output model" contract subsection above.
- **GROMACS: `gmx trjconv` cannot combine `-fit` with `-pbc mol`** — two separate passes, PBC first; this is *why* `fix_PBC.sh` and `strip_and_align_trajectory.sh` are separate scripts.
- **GROMACS: `-fit rot+trans` rotates coordinates but NOT the box** — a fitted frame has a meaningless box, so any PBC-aware selection on it (`gmx select ... within`) silently returns the wrong atoms (measured: 8 spurious waters, one 23 Å out). Order is **PBC → select → fit → write**; the index survives fitting.
- **GROMACS: a solvent-shell selection must come AFTER PBC centering** — `within` is PBC-aware and will otherwise grab periodic images written far from the protein. Sanity-check max water→protein distance ≈ shell + ~1 Å.
- **GROMACS: `gmx select` names a selection with a LEADING QUOTED STRING** — `'"Shell" group "Protein" ...'`, not `'Shell = ...'` (that declares a *variable* → "Too few selections provided"). Dynamic selections also get `_f0_t<time>` appended in the `-on` ndx, so select the group by **index 0**, not by name.
- **Bash: a `[[ … ]] && echo` as the LAST line of a `set -e` script exits 1** — a false test makes SLURM report a successful job as FAILED; use a real `if` (also applies to `((n++))` / `grep -q` as the final statement).
- **GROMACS: never parse a progress readout — `gmx check`'s `Last frame` is throttled AND carriage-return drawn** — on a long trajectory it is never printed at all (2501 frames → stopped at "Reading frame 2000", exit 0), so the parse fails for *some* inputs only. Get the last frame time from the `.gro` title of `trjconv -dump 999999999` instead.
- **Python packaging: setuptools flat-layout auto-discovery fails here** — `gromd_analysis/` is at the repo root next to `scripts/ example/ docs/ knowledgebase/ logs/`, so `pip install -e .` dies with "Multiple top-level packages discovered" unless `pyproject.toml` keeps its explicit `[tool.setuptools] packages = ["gromd_analysis"]`. Install with `--no-deps` (conda owns the scientific deps); `PYTHONPATH=<repo root>` is the no-install fallback.
- **PLUMED/REST2: `partial_tempering` does NOT scale CHARMM CMAP** — measured: `[ cmap ]`/`[ cmaptypes ]` are byte-identical between λ=1.0 and λ=0.5 while charges/LJ/dihedrals scale, so CHARMM+REST2 silently runs a partially-scaled solute Hamiltonian (and both existing self-checks still pass). AMBER has no CMAP and is unaffected; T-REMD/MD scale nothing, so CHARMM is fine there. Same class hits any FF with a term outside the five scaled sections (ff19SB CMAP, Drude polarization) — `docs/FORCE_FIELDS.md` has the per-FF check; run it before REST2 on a new FF.
- **GROMACS: there is no `pdb2gmx -ff list`** — `list` is read as a force-field name and aborts (with a scary `MPI_ABORT`); run `pdb2gmx` with no `-ff` for the menu, or `ls "$GMXLIB"/*.ff`. Confirm a force field resolves with a throwaway `pdb2gmx -ff <name>` run.
- **GROMACS: CHARMM force fields need force-switched vdW (not the AMBER plain cutoff)** — CHARMM36m requires `vdw-modifier=force-switch`, `rvdw-switch=1.0`, `rvdw=rcoulomb=1.2`, `DispCorr=no`; the AMBER cutoff mdp runs but is silently wrong. Selected per job via `FF=charmm36-feb2026_cgenff-5.0` (mdp auto-gated on `FF==charmm*`, forces `CUTOFF_NM=1.2`, AMBER path byte-identical). FF is the MacKerell **force-switch** port at `$HOME/opt/gromacs/ff/charmm36-feb2026_cgenff-5.0.ff` (found via `GMXLIB` in `site_config.sh`; not the `ljpme` port) — **keep the upstream dated directory name**, do not rename it to `charmm36m.ff`, or `parameters.txt` no longer records which release ran. For typing convenience set `FF=charmm36m`: `site_config.sh`'s `FF_ALIASES`/`resolve_ff` expands it at STEP 1 and `parameters.txt` logs the resolved dated name. `WATER=tip3p` auto-resolves to CHARMM-modified TIP3P inside the FF dir. Water always follows the force field — never cross them.
