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

- **The analysis layer is a hybrid: shell drives GROMACS, Python is an installable
  package.** The split line is "does the step invoke `gmx`?" — the ~1200 lines of
  `scripts/analysis/*.sh` stay shell; the ~1000 lines that never touch GROMACS moved
  into `gromd_analysis/` (flat layout at the repo root) with `gromd-*` console
  scripts. *Why not port the shell to Python:* every one of those steps is a `gmx`
  CLI call with a stdin group selection, so a Python port would be `subprocess`
  wrappers around the same commands — one more layer, and it would cost the property
  that each step is echoed as a command you can paste and re-run by hand. *Why
  package the Python:* the parsers were only reachable by `sys.path` hacks, so other
  analyses re-implemented them from scratch (see
  `RELE_simulations/2026-03-25/ai_analysis_test/`, ~1500 lines with hardcoded absolute
  paths). *Why flat and not `src/`:* `PYTHONPATH=<repo root>` then works as a
  no-install fallback; `src/` buys isolation that is worthless for an always-editable
  install. See `GOTCHAS.md` for the auto-discovery trap this creates.

- **`JobDir` (`gromd_analysis/layout.py`) is the one place that knows the output
  layout.** `run_analysis.sh` gets `MODE/TPR/XTC/PREFIX/NCHAINS/...` from
  `eval "$(gromd-layout OUTDIR REP)"` rather than re-deriving them in bash. *Why:*
  the engine/chain-count detection was ~60 lines of shell that every other consumer
  had to duplicate. `JobDir.detect` parses rather than validates — it returns a
  `JobDir` whose paths all exist, or raises `JobLayoutError`; there is no
  half-resolved state to hand downstream.

- **Conformational clustering uses `gromd-cluster` (`gromd_analysis/clustering.py`; sklearn DBSCAN on aligned
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

- **The exchange-statistics block parses into `ExchangeStats`, not a dict.** It was a
  7-key dict that `main()` then mutated with two more keys from a second parser before
  handing it to `report()` 130 lines away. *Why:* the shape is fixed and matters, so it
  belongs in a frozen dataclass with the units and per-pair semantics documented on the
  fields (repo rule: dataclasses, not raw dicts). The two parsers merged into one
  `ExchangeStats.parse`, which also asserts the real invariant the dict version let slide
  — the per-pair lists must have length `n_replicas - 1`, or a truncated block silently
  yields a short table. Parsers raise `RemdLogError` rather than calling `sys.exit`, so
  they are usable as a library; `main()` is the only place that turns that into an exit
  code (same split as `JobLayoutError`). Verified behaviour-preserving: byte-identical
  `remd_acceptance.csv` on all three shipped REMD/REST2 examples.


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
