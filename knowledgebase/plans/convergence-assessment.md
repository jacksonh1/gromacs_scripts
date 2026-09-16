# Plan: convergence assessment for the analysis layer

**Status:** PROPOSED (2026-09-04). Not yet implemented — this document is the design only.

**Scope decided:** all three tools + a rollup, wired into `run_analysis.sh`, with tests.
RMSF split-halves deferred to phase 2 (it is the one piece needing a `gmx` recompute).

## Goal and stance

Give the user a **sense of how converged a trajectory is** — not a yes/no verdict.
Convergence in MD is not provable; only *non*-convergence is detectable. So every tool
here emits **metrics + diagnostic plots** and leaves the judgement to the reader. No tool
prints "converged: true/false".

The assessment characterizes the same designed-structure runs the rest of the analysis layer
does (drift from the designed pose, flexibility, conformational spread), asking the further
question: *given the length of this run, how much can I trust those numbers?*

## Design principle: compose from independent single-concern tools

Three orthogonal questions, three independent tools, plus one rollup. Each tool:

- reads a file **already produced** by `run_analysis.sh` (an `.xvg`, or the clustering
  `assignments.csv`) — so **no new `gmx` invocation** is required,
- is therefore **pure Python** and lives in `gromd_analysis/` (matches the layer's rule:
  shell drives GROMACS, Python does everything else),
- is **layout-blind** and takes explicit input/output paths (the base-script convention),
- installs a `gromd-*` console entry point in `pyproject.toml`,
- is runnable standalone and re-runnable by hand.

This mirrors how the existing metrics are split — one script per metric — so a user can run
just the piece they care about.

## The three concerns

### 1. Scalar-timeseries convergence — `gromd-ts-convergence`
New module `gromd_analysis/ts_convergence.py`.

**Input:** one scalar `.xvg` (single y-column vs time), parsed with the existing
`parse_xvg`. Run on `<prefix>_rmsd.xvg` and `<prefix>_rg.xvg`.

**Computes (per series):**
- **Equilibration onset `t0`** — automatic, choose `t0` that maximizes the effective sample
  size of the remaining data (reverse-cumulative averaging, Chodera 2016). Reports how much
  of the front of the run is discarded as equilibration. All stats below use `t > t0`.
- **Integrated autocorrelation time `τ_int`** and **effective sample size
  `N_eff = N / (2·τ_int)`**. This is the headline number: how many *independent* samples the
  run actually contains, which is what limits every mean/error, not the raw frame count.
- **Block-averaged SEM** (Flyvbjerg–Petersen): standard error vs block size; report the
  plateau value and the block size at which it plateaus. Cross-checks `τ_int`.
- **Trailing-window drift** — linear slope of the observable over the last window, with a
  confidence interval. A CI that excludes zero means the observable is still drifting. This
  reuses the plateau-slope *idea* from `scripts/simulation/density_converged.py` (kept
  consistent with it; not imported — that lives in the simulation layer and is stdlib-only).

**Output:** `<prefix>_<metric>_convergence.json` (all numbers) and
`<prefix>_<metric>_convergence.png` (the series with the `t0` marker and running mean;
a block-SEM-vs-block-size inset; the trailing slope ± CI).

### 2. Distribution stationarity — `gromd-split-halves`
New module `gromd_analysis/split_halves.py`.

**Input:** the same scalar `.xvg` (RMSD / Rg).

**Computes:** compares the first half vs the second half of the post-`t0` data —
- KS statistic and histogram overlap coefficient,
- per-half mean / SD, and the mean difference expressed in units of the block-SEM.

Two halves whose distributions disagree is a direct, intuitive sign the run has not
converged. Cheap complement to `τ_int`.

**Output:** `<prefix>_<metric>_halves.json` and `_halves.png` (overlaid half histograms +
a stats box).

**RMSF split-halves — phase 2.** A per-residue "did flexibility change over the run" check
needs RMSF recomputed on each half (`gmx rmsf -b/-e`), i.e. a `gmx` call, so it does **not**
belong in this pure-Python tool. Add later as a shell `calc_rmsf_halves.sh` (two `gmx rmsf`
passes) plus a Python comparison (per-residue Pearson r and |Δ RMSF|), which also flags
*which* residues destabilize late in the run.

### 3. Conformational-space exploration — `gromd-sampling-convergence`
New module `gromd_analysis/sampling_convergence.py`.

**Input:** the clustering `clustering/<prefix>_cluster_assignments.csv` already written by
`gromd-cluster`. No `gmx`.

**Computes:**
- **State-discovery curve** — cumulative count of distinct clusters vs time. If it is still
  climbing at the end of the trajectory, conformational space is still being explored.
- **Population stationarity** — cluster populations in the first vs second half, compared
  with the Jensen–Shannon divergence.
- **Noise fraction** — fraction of frames DBSCAN labelled `-1`. A high value means the run
  is under-sampled / heterogeneous relative to the clustering cutoff.

**Output:** `<prefix>_sampling_convergence.json` and `_sampling_convergence.png`.

## Rollup

`<prefix>_convergence_summary.txt` — one table of the headline numbers across all three
tools (`t0`, `τ_int`, `N_eff`, drift slope ± CI, half-split KS, whether the state-discovery
curve has plateaued, noise %), followed by plain-language guidance on how to read them.
**No verdict line.**

Output format decision: **JSON per tool** (machine-readable, importable from a notebook via
the package, consistent with library use) **plus the human `summary.txt` rollup**.

## Integration — `run_analysis.sh`

Append `[CMD]` steps after the producers each tool consumes:
- `ts-convergence` and `split-halves` after the RMSD and Rg `.xvg` steps,
- `sampling-convergence` after the clustering step.

All non-fatal (`|| echo "[WARN] …"`), exactly like every other step, so a convergence-tool
failure never fails the job. Same path for both pipelines; for REMD the tools run on the
**rep000** (lowest-T, 300 K constant-temperature) outputs, identical to the other metrics.

## REMD caveats (also add to `GOTCHAS.md` + `scripts/analysis/README.md`)

- Run the assessment on **rep000** — the 300 K constant-temperature slot. Higher slots are
  higher-temperature ensembles and are not the quantity of interest.
- For REMD, `τ_int` measured on rep000 is legitimately **short**: exchanges decorrelate the
  slot every `REPLEX_PS`, so `N_eff` from rep000 is real and often large. But a converged
  REMD result *also* requires good replica mixing, which is a different question already
  answered by `gromd-acceptance` (per-pair exchange rates) and `gromd-roundtrip` (round
  trips / dwell). The new tools **complement** those; they do not replace them. A run with
  poor mixing can still show a deceptively short `τ_int` in rep000.

## Tests (`tests/`)

The metrics are pure Python, so they are unit-testable against synthesized series in
`tmp_path` (never reading `example/outputs/`, per the project rule):
- a stationary series vs a linearly drifting one → assert the drift series' trailing-slope CI
  excludes zero and its `N_eff` is lower;
- a strongly autocorrelated series → assert `N_eff` ≪ N;
- a bimodal (shifted second half) series vs a unimodal one → assert the split-halves KS
  statistic separates them;
- a synthetic `assignments.csv` where new clusters appear only late → assert the
  state-discovery curve has not plateaued.

## Files this will touch (when built)

- new: `gromd_analysis/ts_convergence.py`, `gromd_analysis/split_halves.py`,
  `gromd_analysis/sampling_convergence.py`
- new: `gromd_analysis/convergence_summary.py` (rollup) — or fold the rollup into
  `run_analysis.sh` if it is only assembling the three JSONs
- edit: `pyproject.toml` (four `gromd-*` entry points), `scripts/analysis/run_analysis.sh`
  (the new `[CMD]` steps), `scripts/analysis/README.md` (document the tools + REMD caveat),
  `knowledgebase/GOTCHAS.md` + repo `CLAUDE.md` (REMD rep000 / `τ_int` caveat),
  `tests/` (one test module per tool)

## Dependencies

numpy (already present); **scipy** for the KS test and Jensen–Shannon divergence. If scipy is
not already in `environment.yml`, **add it** — do not hand-roll these standard statistics (see
the no-hand-rolling rule in `09-fragfold/CLAUDE.md`). scipy is a mainstream, well-tested
dependency; a bespoke KS / JS implementation is an unnecessary source of subtle numerical bugs.
