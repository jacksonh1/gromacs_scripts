# Analysis Scripts

GROMACS post-processing utilities, shared by the T-REMD and plain-MD pipelines.

These metrics characterize how a designed structure behaves over the trajectory:
**RMSD** vs the reference = drift from the designed pose (stability / design retention,
objectives #1 and #3); **RMSF** = per-residue flexibility (flexible-region identification,
#2); **Rg** and **DSSP** track overall compactness and secondary-structure persistence.
All are measured relative to the starting (designed) structure — runs never begin from an
unfolded state.

## Design: layout-blind scripts

The base scripts take **explicit input/output file paths** and know nothing
about directory layout — the caller (each sbatch, or you) supplies the paths.
This is why one set of scripts serves both pipelines: the only difference
between them is *where files live* (`prod/rep000/remd.*` vs `prod/md.*`).

The single exception is `fix_PBC_strip_align.sh`, the orchestrator: it is the
one place that encodes the `<prefix>_pbc.xtc` intermediate naming convention.

## Design: shell drives GROMACS, Python does the rest

The analysis layer is split along one line — **does the step invoke `gmx`?**

- **Shell (`scripts/analysis/*.sh`)** — every step that drives the GROMACS CLI:
  PBC correction, stripping, alignment, `gmx rms` / `gyrate` / `rmsf` / `dssp`,
  solvated-snapshot extraction. These stay as shell because they *are* `gmx`
  invocations with stdin group selections; wrapping them in `subprocess` would
  add a layer and cost the property that every step is echoed as a command you
  can paste and re-run by hand.
- **Python (`gromd_analysis/`, at the repo root)** — everything that never touches
  GROMACS: parsers, clustering, plotting, job-directory detection. This half is an
  installable package, so it is importable from notebooks and other repos instead of
  being copy-pasted or re-implemented.

| command | module | does |
|---|---|---|
| `gromd-layout` | `layout.py` | parse OUTDIR → `JobDir` (mode, tpr, xtc, prefix, chain count) |
| `gromd-plot-xvg` | `xvg.py` | parse + plot any GROMACS `.xvg` |
| `gromd-plot-dssp` | `dssp.py` | secondary-structure map |
| `gromd-cluster` | `clustering.py` | conformational clustering (sklearn) |
| `gromd-acceptance` | `remd_log.py` | REMD/REST2 exchange acceptance rates |
| `gromd-chain-index` | `chains.py` | per-chain `.ndx` groups from the topology |

Install (once, into the analysis env — `install_python_env.sh` does this for you):

```bash
conda activate groMD_env && pip install --no-deps -e .   # from the repo root
```

`--no-deps` because conda already provides matplotlib/numpy/MDAnalysis/scikit-learn
from `environment.yml`; pip must not reinstall and shadow them. On a node without
network, `PYTHONPATH=<repo root>` makes `import gromd_analysis` work with no install
(the package uses a flat layout for exactly this reason) — but the `gromd-*` commands
need the install.

### Library use

`JobDir` is the entry point for any one-off analysis, so a notebook never has to
hardcode `prod/rep000/remd.xtc` or re-derive the chain count:

```python
from pathlib import Path
from gromd_analysis import JobDir, parse_xvg

job = JobDir.detect(Path("example/outputs/output_T-REMD/<run>"))
job.mode, job.n_chains, job.xtc          # ('REMD', 1, PosixPath('.../prod/rep000/remd.xtc'))
rmsd = parse_xvg(Path(f"{job.prefix}_rmsd.xvg"))   # prefix is a stem, not a file
```

`JobDir.detect` either returns a `JobDir` whose paths all exist, or raises
`JobLayoutError` — it never hands back a half-resolved directory.

### Tests

```bash
conda activate groMD_env && pytest        # from the repo root
```

`tests/` covers the parsers — `.xvg` (including the nm→Å relabelling), the topology
`[ molecules ]` / `[ moleculetype ]` sections, `JobDir.detect` for all three engines
plus its failure modes, the replica-exchange statistics block, and the DSSP reader.
Every fixture is synthesized in `tmp_path`; **no test reads `example/outputs/`**,
which is `.gitignore`d and whose bulk stage dirs are scratch symlinks that get purged.

Only the GROMACS-free code is covered. The `.sh` steps are verified by running
`run_analysis.sh` on a finished job and diffing the outputs — the shipped examples
are there for exactly that.

---

## Trajectory preparation

### `fix_PBC.sh` — PBC correction (full system)
```
bash fix_PBC.sh TPR XTC OUT_PBC_XTC
```
Fixes periodic boundary condition artifacts using `gmx trjconv -pbc mol -center -ur compact`.
Output: a full-system (protein + water + ions) PBC-corrected trajectory at `OUT_PBC_XTC`.

This is the **single-chain** path. For a multi-chain complex `-pbc mol` can split
the chains across a box boundary — see "Multi-chain complexes" below.

Use this directly when you need to inspect the trajectory with solvent included.
For the typical automated case (protein analysis only), use `fix_PBC_strip_align.sh`.

**`-pbc nojump` is intentionally not used.** `nojump` compares consecutive
frames to detect box-crossing events, but T-REMD coordinate exchanges cause
discontinuous jumps between frames that `nojump` misinterprets. `-pbc mol` is a
per-frame operation and is safe for REMD trajectories.

---

### `strip_and_align_trajectory.sh` — strip waters/ions + backbone alignment
```
bash strip_and_align_trajectory.sh TPR PBC_XTC OUT_PREFIX [REF_GRO]
```
Requires a PBC-corrected, full-system trajectory (run `fix_PBC.sh` first). Fits
each frame to a reference (default: first frame) using backbone atoms (N, CA, C,
O) and outputs only protein atoms. Writes:
- `<OUT_PREFIX>_stripped_aligned.xtc` — protein-only, backbone-aligned trajectory (kept)
- `<OUT_PREFIX>_stripped_aligned.gro` — protein-only reference structure (kept; the reference
  for RMSD, clustering, etc.)
- `<OUT_PREFIX>_frame0_fullsys.gro` — full-system frame 0, used internally as `-s` to
  match the XTC atom count, then deleted

The `-s` reference must match the XTC atom count (see "atom-count matching" in
`REMD_log_reference.md`); the throwaway full-system GRO exists for exactly that.

---

### `fix_PBC_strip_align.sh` — PBC fix + strip + align, no intermediate file
```
bash fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]
```
Orchestrator: runs `fix_PBC.sh` (writing `<OUT_PREFIX>_pbc.xtc`) then
`strip_and_align_trajectory.sh`, then deletes the large `_pbc.xtc` intermediate.
Same kept outputs as the two-script sequence: `<OUT_PREFIX>_stripped_aligned.xtc` and
`<OUT_PREFIX>_stripped_aligned.gro`.

Use this in the automated pipeline. Use `fix_PBC.sh` directly when you need to
keep the solvent-included trajectory.

---

### `extract_protein.sh` — protein-only structure from a full-system frame
```
bash extract_protein.sh STRUCT TPR OUT_GRO
```
Selects the `Protein` group from a full-system structure (e.g. `em/em.gro`) and
writes it on its own. The output has the same atom ordering as the protein-only
trajectory from `strip_and_align_trajectory.sh`, so it works directly as a
reference for `calc_traj_rmsd.sh`.

Typical use: build an **initial-structure reference** from the minimized `em.gro`
so RMSD measures drift from the starting structure rather than from the first
frame of the trajectory. Both pipelines do this automatically (writing
`<prefix>_init.gro`) and pass it as the RMSD reference.

The extraction runs `trjconv -pbc whole`, so the protein is made whole even if it
straddles a box boundary in `em.gro`. This matters because `gmx rms` does **not**
make its reference whole — a broken reference silently corrupts every RMSD value
(in one example it inflated the t=0 backbone RMSD to ~13 Å; with the whole
reference it reads a physical ~3 Å). For a multi-chain complex the reference also
needs the chains kept in one image — see "Multi-chain complexes" below.

---

## Solvated snapshots: `extract_solvated_snapshots.sh` (MD only)

```
bash extract_solvated_snapshots.sh TPR XTC OUT_PREFIX [SHELL_NM] [N_SNAPSHOTS]
```

Writes a **handful of PBC-corrected PDBs that keep a solvent shell** — the protein
plus every whole water/ion within `SHELL_NM` (default `0.5` nm = 5 Å) of it. Defaults
to 5 snapshots evenly spaced over the trajectory (first, three middle, last), all
fitted onto a common reference so they superimpose when loaded together.

Everything else in this directory works on the *protein-only* trajectory; this is the
one output that keeps solvent, for looking at interface/bound-state waters.

Outputs go to a **top-level `solvated_snapshots/` dir in `OUTDIR`** (a sibling of
`analysis/`, not inside it):

```
<OUTDIR>/solvated_snapshots/md_solvshell_000000ps.pdb
<OUTDIR>/solvated_snapshots/md_solvshell_000500ps.pdb   ... etc
```

Zero-padded ps, so lexical order is chronological.

**Why standalone PDBs and not a trajectory.** A shell selection is *dynamic*: the set
of waters within `SHELL_NM` changes every frame, so the atom count changes every frame.
XTC/TRR require a fixed atom count and no single topology would match. Independent
one-frame PDBs sidestep that — each carries its own atoms. This is also why the count
is deliberately small: a few files you can open at once, not hundreds.

**Order of operations — all three constraints are real:**

1. **PBC before selection.** `gmx select`'s `within` is PBC-aware and will select a
   water whose *periodic image* is near the protein while its stored coordinate sits
   across the box — written out as a water floating in space.
2. **Selection before fitting.** `-fit rot+trans` rotates coordinates but leaves the
   **box vectors untouched**, so on a fitted frame the PBC-aware `within` measures
   against a box that no longer matches. Measured on the 2 ns example: 1622 atoms
   selected on the fitted frame vs 1598 pre-fit — 8 bogus waters, some 23 Å out.
   Atom numbering is unchanged by fitting, so the index applies to the fitted frame.
3. **PBC and fitting must be separate `trjconv` calls** — `trjconv` refuses `-fit`
   together with `-pbc mol`. (Same reason `fix_PBC.sh` and
   `strip_and_align_trajectory.sh` are separate scripts.)

`run_analysis.sh` runs this automatically for **MD jobs only** and dispatches on chain
count like everything else. Override the defaults via env vars:
`SHELL_NM=0.8 N_SNAPSHOTS=9 bash run_analysis.sh OUTDIR`.

To add snapshots to an *older* job that finished before this step existed, call the
script directly rather than re-running `run_analysis.sh` — the latter also redoes the
full PBC/strip/align pass over the whole trajectory and overwrites the existing
`analysis/` outputs, while this reads only the frames it dumps.

---

## Trajectory metrics (trajectory-agnostic)

All four share one signature and operate on the protein-only, aligned outputs of
`fix_PBC_strip_align.sh`:
```
bash calc_traj_<metric>.sh STRUCT TRAJ OUT
#   STRUCT = protein-only reference .gro  (e.g. analysis/md_stripped_aligned.gro)
#   TRAJ   = protein-only aligned .xtc    (e.g. analysis/md_stripped_aligned.xtc)
#   OUT    = output file (.xvg, or .dat for dssp)
```
STRUCT and TRAJ are both protein-only, so their atom counts match and STRUCT is
safe as `-s`.

| Script | Computes | Output |
|--------|----------|--------|
| `calc_traj_rmsd.sh` | Backbone RMSD vs time (`gmx rms`) | `.xvg` |
| `calc_traj_rg.sh`   | Radius of gyration vs time (`gmx gyrate`) | `.xvg` |
| `calc_traj_rmsf.sh` | Per-residue backbone RMSF (`gmx rmsf -res`) | `.xvg` |
| `calc_traj_dssp.sh` | DSSP secondary structure over time (`gmx dssp`) | `.dat` |

Each QC/metric is a separate script so they can be modified independently.

**RMSD reference (`STRUCT` for `calc_traj_rmsd.sh`).** `gmx rms` re-fits every frame
to `STRUCT` before measuring, so (a) `STRUCT` *is* the RMSD reference, and (b) the
reference does **not** need to be pre-aligned to the trajectory — the fit is internal
and per-frame. The pipelines pass `<prefix>_init.gro` (the minimized initial structure,
from `extract_protein.sh`, made whole) so RMSD reports drift from the starting structure.
To instead measure frame-to-frame drift during the run, pass `<prefix>_stripped_aligned.gro`
(the first frame). The other metrics use `_stripped_aligned.gro` as topology only.

---

## Plotting (standalone, reusable)

Decoupled from computation — they render data files already on disk, so you can
re-plot/restyle without re-running `gmx`. Require `matplotlib`.

### `gromd-plot-xvg` — generic line plot
```
gromd-plot-xvg PLOT_XVG OUT_PNG
```
Reads the title and axis labels embedded in the `.xvg` (`@ title` / `@ xaxis
label` / `@ yaxis label`) and plots every data column against column 0. Serves
the RMSD, Rg, and RMSF outputs — and any other xvg — unchanged.

GROMACS writes distances in nm; any axis labelled in `nm` is converted to **ångström**
(×10, label rewritten to `Å`) so the plots read in the units expected for structural
work. The underlying `.xvg` data files are left in nm.

### `gromd-plot-dssp` — secondary-structure map
```
gromd-plot-dssp DSSP_DAT OUT_PNG
```
Renders the `gmx dssp` `.dat` as a residue × frame categorical heatmap with a
legend. Colours are a fixed semantic palette (helices in blues/purple, strands in
warm reds, turns/bends in greens, coil in light grey) so the map is consistent
across runs and structured regions stand out against the loopy background.

---

## Conformational clustering: `gromd-cluster`
```
gromd-cluster STRUCT TRAJ OUT_PREFIX [--method dbscan|kmeans]
    [--cutoff NM] [--min-samples N] [--n-clusters K] [--selection SEL] [--stride N]
```
Groups the sampled frames into discrete conformational states (Cα clustering) and
writes, per state, its population and a representative structure. Operates on the same
protein-only, PBC-fixed, aligned trajectory as the metrics above, so it serves both
pipelines and runs once in the shared part of `run_analysis.sh` (no multi-chain variant).

**Why not `gmx cluster`.** The gromos algorithm builds the full pairwise-RMSD matrix —
**O(N²)** in time and memory, impractical for the long production runs (25k+ frames).
`gromd-cluster` clusters on flattened Cα coordinates with scikit-learn, which scales:
a 25k-frame run clusters in ~13 s. `-pbc cluster` in the multi-chain pipeline is a
*periodic-image* operation and is **unrelated** to this conformational clustering.

**Methods.**
- **DBSCAN** (default): density-based; finds the cluster count automatically and labels
  sparsely-visited frames as *noise* (`-1`). The `--cutoff` is a backbone-RMSD cutoff in
  **nm** (default `0.20` = 2.0 Å). Because the input frames are pre-aligned to one common
  reference, the flattened-coordinate Euclidean distance equals `√N_atoms × RMSD`, so the
  script sets the DBSCAN radius `eps = cutoff_Å × √N_selected` — i.e. `--cutoff` is a true
  RMSD threshold. (This is the one fix over the Amber `cluster_MD.py`, whose `eps` was a
  raw flattened-coord distance mislabelled as Å.)
- **k-means** (`--method kmeans --n-clusters K`): O(N), always scales; forces every frame
  into one of K clusters with no noise.

**Cluster count vs noise — `--min-samples`.** This is the DBSCAN density knob: how many
frames must lie within the cutoff of a core point for a region to count as a cluster
rather than noise. A fixed value doesn't scale (10 frames is 5% of a 200-frame run but
0.1% of a 10k-frame run → a long tail of sub-1% clusters), so the **default is adaptive:
`max(10, 1.5% of frames)`** — a region must hold ~1.5% of the trajectory to be a state;
sparser frames fall into noise. This keeps even a heterogeneous run to ≲10 clusters (tuned
on the WW-domain REMD slots). Raise `--min-samples` for fewer/denser clusters; pass an
absolute integer to override the adaptive default (e.g. `--min-samples 10` restores the
old, fragmented behaviour). The effective value is printed and recorded in the summary.

Clusters are relabelled by descending population, so `c00` is the dominant state.
**Outputs** go in a **`clustering/` subdir** next to the prefix (a run can produce many
rep PDBs), named `<prefix-name>_cluster_*`:
e.g. `analysis/remd_rep000` → `analysis/clustering/remd_rep000_cluster_{assignments.csv,
rep_c00.pdb, rep_c01.pdb, …, populations.png, timeseries.png, summary.txt}`.

**Scaling caveat.** DBSCAN's neighbour graph grows large when nearly all frames fall
within the cutoff (a very stable structure at a loose cutoff — the trivial "one cluster"
answer); it still completes but uses more memory at tens of thousands of frames.
`--stride` and `--method kmeans` are the escape hatches.

`run_analysis.sh` runs DBSCAN at `CLUSTER_CUTOFF` (env override, default 0.20 nm); the
adaptive `min_samples` applies automatically.

---

## `gromd-acceptance` — exchange acceptance rates (REMD only)
```
gromd-acceptance OUTDIR [--rep REP] [--plot]
```
Parses the "Replica exchange statistics" block from the GROMACS log (pre-computed
by GROMACS at the end of each run) and reports per-pair empirical acceptance
rates, mean Metropolis probabilities, and exchange counts. Outputs a console
table and `OUTDIR/analysis/remd_acceptance.csv`. Use `--plot` to also write
`remd_acceptance.png`.

Target acceptance rates for T-REMD are 20–30% per adjacent pair. Rates above 50%
indicate replicas are too closely spaced in temperature.

---

## Re-running the whole analysis: `run_analysis.sh`

```
bash run_analysis.sh OUTDIR [REP]
```
One entry point for the entire post-analysis. Auto-detects MD vs T-REMD from the
job layout (`prod/md.tpr` vs `prod/rep<REP>/remd.tpr`), activates the Python env,
runs the full chain (acceptance rates for REMD → PBC fix/strip/align → whole RMSD
reference → RMSD/Rg/RMSF/DSSP + plots → conformational clustering), and echoes each
`[CMD]` as it goes. `REP`
defaults to `000` and is ignored for plain MD.

This is exactly what each engine's sbatch calls as its post-analysis step, so it
is the canonical way to **re-run analysis without resubmitting the simulation** —
e.g. after editing an analysis script:
```bash
bash run_analysis.sh OUTDIR          # plain MD
bash run_analysis.sh OUTDIR 000      # T-REMD, lowest-T slot
```
Re-runs overwrite cleanly (GROMACS backups are disabled inside the script), so it
is idempotent. It needs the production trajectory still present — under the
folder-symlink output model it reads `prod/md.xtc` or `prod/rep<REP>/remd.xtc`, which
resolves through the `prod/` folder symlink into scratch (or is a real file when
`SYMLINK_BULK=0`). Copy anything you need off scratch before it is purged. Jobs
predating the output restructure kept their trajectory in a top-level
`trajectories/` dir; that layout is accepted (with an `[INFO]` line saying so) purely
so old job dirs stay analysable — nothing writes there any more.

Knobs are environment variables, not job parameters, because this script is meant to
be re-run by hand — that is where you would want to vary them:

| Variable | Default | Effect |
|---|---|---|
| `CLUSTER_CUTOFF` | `0.20` | backbone-RMSD clustering cutoff, nm |
| `SHELL_NM` | `0.5` | solvated-snapshot solvent shell, nm (MD only) |
| `N_SNAPSHOTS` | `5` | number of solvated snapshots (MD only) |

```bash
SHELL_NM=0.8 N_SNAPSHOTS=10 bash run_analysis.sh OUTDIR
```

### On a compute node: `analysis.sbatch`

`run_analysis.sh` is fine on a login node for a small job, but a long trajectory
(GB-scale XTC, clustering over thousands of frames) belongs in a SLURM allocation:

```
sbatch --export=ALL,GROMACS_SCRIPTS_DIR=<repo>/scripts/simulation,OUTDIR=<job> analysis.sbatch
```
The usual entry point is to copy **`example/submit_jobs/submit_analysis.sh`**, set
`OUTDIR`, and run it — same copy-and-edit pattern as `submit_MD.sh` / `submit_REMD.sh`.

`analysis.sbatch` supplies only the allocation and the module/GROMACS environment;
all the work is still `run_analysis.sh`, so MD/T-REMD/REST2 detection and the knobs
above behave identically. It requests **no GPU** — every analysis step is CPU-only,
so on `pi_keating` it lands on the CPU-only nodes instead of queueing behind
production runs.

## Multi-chain complexes

`run_analysis.sh` counts the protein chains in the topology (`build/*.top`,
`[ molecules ]`) and **dispatches automatically**:

- **1 chain** → the standard scripts above, unchanged.
- **>1 chain** → a separate `multichain_*` pipeline. Single-chain runs never touch it.

Why a separate path: for a complex, `-pbc mol` wraps each chain's COM independently
and `-pbc whole` only un-breaks within a molecule — so the chains can land in
different periodic images and the complex is **split** across a box boundary,
corrupting the RMSD reference (observed: a split reference gave a median backbone
RMSD of 24.5 Å while Rg stayed ~1.4 nm). The multi-chain scripts apply
`gmx trjconv -pbc cluster`, which pulls the chains into **one image**. It is
per-frame (REMD-safe) and, for a single chain, byte-identical to `-pbc whole`.

> **Naming:** these scripts are `multichain_*` and use `-pbc cluster` only as a
> *periodic-image* fix. That is unrelated to **conformational clustering**, which lives
> in `gromd-cluster` (sklearn DBSCAN, see above) — not `gmx cluster`. The PBC scripts
> are deliberately *not* named `*cluster*` to keep that distinction clear.

The multi-chain scripts (mirror their single-chain counterparts):

| Script | Mirrors | Role |
|--------|---------|------|
| `multichain_fix_PBC.sh` | `fix_PBC.sh` | PBC fix in 3 passes: `whole → cluster → mol+center+compact` |
| `multichain_fix_PBC_strip_align.sh` | `fix_PBC_strip_align.sh` | orchestrator (reuses `strip_and_align_trajectory.sh`) |
| `multichain_extract_protein.sh` | `extract_protein.sh` | RMSD reference via `-pbc cluster` |
| `multichain_extract_solvated_snapshots.sh` | `extract_solvated_snapshots.sh` | solvated snapshot PDBs, chains kept in one image |
| `gromd-chain-index` | — | chain detector + index builder (per-chain + per-chain-backbone groups) |
| `multichain_chain_rmsd.sh` | `calc_traj_rmsd.sh` | one chain's backbone RMSD, fit to itself |
| `multichain_chain_rmsf.sh` | `calc_traj_rmsf.sh` | one chain's per-residue RMSF |
| `multichain_interchain_dist.sh` | — | inter-chain minimum-image distance (`gmx mindist`) |

**Output mirrors the single-chain set** (`<p>_stripped_aligned.{xtc,gro}`,
`<p>_init.gro`, `<p>_{rmsd,rg,rmsf,dssp}.*` on the whole complex), **plus**:
`<p>_chain{A,B,…}_{rmsd,rmsf}.{xvg,png}` and `<p>_interchain_mindist.{xvg,png}`
(one curve per chain pair). Per-chain RMSD fits each chain to *itself*, so it reports
that chain's internal drift regardless of how the chains sit relative to each other.

`gromd-chain-index` uses the protein **molecule types** from the topology
(the physically-correct chains), not `gmx splitch` (which over-split one observed
system 2→3 on an internal residue-numbering gap).

**Bound-complex assumption.** `-pbc cluster` (and whole-complex RMSD) assume the
chains stay within ~half the box. The `interchain_mindist` curve is the check: if it
climbs toward half the box, the complex is dissociating *and* the box is too small
(minimum-image violation) — a setup problem, not just an analysis one.

## Typical workflows (manual / per-metric)

`run_analysis.sh` is the normal path; the calls below are the individual steps it
runs, for when you want just one metric or a custom reference.

### Plain MD
```bash
TPR=OUTDIR/prod/md.tpr
XTC=OUTDIR/prod/md.xtc
P=OUTDIR/analysis/md

bash   fix_PBC_strip_align.sh "$TPR" "$XTC" "$P"
bash   calc_traj_rmsd.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rmsd.xvg" && gromd-plot-xvg "${P}_rmsd.xvg" "${P}_rmsd.png"
bash   calc_traj_rg.sh   "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rg.xvg"   && gromd-plot-xvg "${P}_rg.xvg"   "${P}_rg.png"
bash   calc_traj_rmsf.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rmsf.xvg" && gromd-plot-xvg "${P}_rmsf.xvg" "${P}_rmsf.png"
bash   calc_traj_dssp.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_dssp.dat" && gromd-plot-dssp "${P}_dssp.dat" "${P}_dssp.png"
```

### T-REMD
Same as above, but the inputs are the lowest-T slot and the prefix is
`remd_rep000`, plus the REMD-only acceptance-rate QC:
```bash
TPR=OUTDIR/prod/rep000/remd.tpr
XTC=OUTDIR/prod/rep000/remd.xtc
P=OUTDIR/analysis/remd_rep000

gromd-acceptance OUTDIR                       # QC: exchange acceptance rates
bash   fix_PBC_strip_align.sh "$TPR" "$XTC" "$P"
# ... then the same calc_traj_*/plot_* calls as above
```

---

## Key difference from Amber T-REMD

In GROMACS T-REMD, each replica runs at a **fixed temperature** and **coordinates**
are exchanged between neighboring replicas. This means:

- `prod/rep000/remd.xtc` **is** the lowest-T constant-temperature trajectory — use it directly.
- `prod/rep001/remd.xtc`, etc. are the higher-temperature trajectories.
- **No demux step is needed** to obtain the constant-temperature ensemble.

`gmx demux` generates a trajectory that follows a specific *configuration* as it
walks through temperature space — useful for visualizing the random walk of a
molecule across replicas, but not the correct tool for obtaining the
thermodynamic ensemble at a given temperature.

---

Planned tools:
- Empirical transition matrix analysis (mixing quality between temperature slots)
- Round-trip counter (requires per-frame walker trajectory reconstruction from `Repl ex` lines)
