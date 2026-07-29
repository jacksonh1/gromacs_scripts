# Decisions

The *why* behind non-obvious design choices. One entry per decision. Record the
alternative that was rejected and why. Operational pitfalls live in `CLAUDE.md`
("Known gotchas"), not here — this is rationale, not runbook.

---

## Architecture

- **One engine per method, not branches.** `REMD-gromacs.sbatch`,
  `MD-gromacs.sbatch`, and (in progress) `REST2-gromacs.sbatch` are separate
  scripts. New case-specific capability goes in a parallel pipeline with a
  dispatch, never as `if` branches inside a working engine. *Why:* isolate
  complexity from the paths that already work.

- **Multi-chain is its own `multichain_*` analysis pipeline**, selected
  automatically by `run_analysis.sh` when the topology has >1 protein chain. It
  uses the `-pbc cluster` periodic-image fix (REMD-safe, per-frame). *Why:* keep
  the single-chain path simple; a broken complex image silently inflates RMSD/Rg.
  (Scripts are named `multichain_*`, never `*cluster*`, to avoid confusion with
  conformational clustering.)

- **Stage folders named by role, not ensemble** (`build/ em/ heat/ density/
  relax/ prod/` for MD; `build/ em/ density/ equil/ prod/` for REMD). *Why:* the
  name stays correct regardless of which ensemble runs in it. Do not reintroduce
  `nvt/ npt/ eq/`.

- **Production basenames are method-named** (`md.*`, `remd.*`, `rest2.*`), not
  ensemble-named. *Why:* the analysis layer keys off the basename to auto-detect
  the engine.

## Analysis

- **Conformational clustering uses `cluster_traj.py` (sklearn DBSCAN on aligned
  Cα), not `gmx cluster`.** *Why:* gmx cluster builds the full O(N²) pairwise-RMSD
  matrix — impractical for 25k+ frame production runs. DBSCAN on flattened aligned
  coords scales (~13 s for 25k frames). See `[[reference_clustering]]` /
  `CLAUDE.md`.

- **Constant-temperature ensembles come straight from the slot trajectory** —
  `prod/rep000/remd.xtc` *is* the 300 K ensemble. *Why:* each `rep{i}/` is a fixed
  temperature slot; coordinates exchange between slots. No `gmx demux` needed
  (demux follows a configuration's random walk, not the thermodynamic ensemble).

- **REMD acceptance rates are parsed from the log's pre-computed statistics
  block**, not by recounting per-frame `Repl ex` lines. *Why:* the block already
  holds per-pair acceptance rates + counts.

## Robustness

- **Parameters validated fail-loud at STEP 1** (`validate_params.py`): unknown
  config keys rejected as typos, required/type/range checked on resolved values.
  *Why:* a bad parameter should exit with `[ERROR]` — not silently use a default
  and waste compute. Three lists must stay in sync when adding a param (see
  `CLAUDE.md`).

- **`gmx_mpi` everywhere** (this build ships no plain `gmx`/`gmx mdrun`). Serial
  prep steps run as a single MPI rank without `mpirun`.

## REST2 (in progress — see `plans/REST2-pipeline.md`)

- **Used exactly like REMD: single-chain, whole-protein hot region.** No
  complexes, no residue-subset selection. *Why:* matches how the user runs it;
  keeps `mark_hot_region.py` trivial (mark the whole `Protein*` moleculetype).

- **`rest2.*` basename + a third `run_analysis.sh` mode**, rather than reusing
  `remd.*`. *Why:* rep000 (λ=1) analysis is identical to REMD's, but a
  self-describing basename beats a folder that lies about which method ran. Chosen
  over the zero-code-change option of reusing `remd.*`.
