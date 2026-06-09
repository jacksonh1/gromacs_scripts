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

### `plot_xvg.py` — generic line plot
```
python plot_xvg.py PLOT_XVG OUT_PNG
```
Reads the title and axis labels embedded in the `.xvg` (`@ title` / `@ xaxis
label` / `@ yaxis label`) and plots every data column against column 0. Serves
the RMSD, Rg, and RMSF outputs — and any other xvg — unchanged.

GROMACS writes distances in nm; any axis labelled in `nm` is converted to **ångström**
(×10, label rewritten to `Å`) so the plots read in the units expected for structural
work. The underlying `.xvg` data files are left in nm.

### `plot_dssp.py` — secondary-structure map
```
python plot_dssp.py DSSP_DAT OUT_PNG
```
Renders the `gmx dssp` `.dat` as a residue × frame categorical heatmap with a
legend. Colours are a fixed semantic palette (helices in blues/purple, strands in
warm reds, turns/bends in greens, coil in light grey) so the map is consistent
across runs and structured regions stand out against the loopy background.

---

## Conformational clustering: `cluster_traj.py`
```
python cluster_traj.py STRUCT TRAJ OUT_PREFIX [--method dbscan|kmeans]
    [--cutoff NM] [--min-samples N] [--n-clusters K] [--selection SEL] [--stride N]
```
Groups the sampled frames into discrete conformational states (Cα clustering) and
writes, per state, its population and a representative structure. Operates on the same
protein-only, PBC-fixed, aligned trajectory as the metrics above, so it serves both
pipelines and runs once in the shared part of `run_analysis.sh` (no multi-chain variant).

**Why not `gmx cluster`.** The gromos algorithm builds the full pairwise-RMSD matrix —
**O(N²)** in time and memory, impractical for the long production runs (25k+ frames).
`cluster_traj.py` clusters on flattened Cα coordinates with scikit-learn, which scales:
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

## `remd_acceptance.py` — exchange acceptance rates (REMD only)
```
python remd_acceptance.py OUTDIR [--rep REP] [--plot]
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
is idempotent. It needs the production trajectory still present — `trajectories/*.xtc`
are symlinks into scratch, so copy them off scratch before it is purged.

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
> in `cluster_traj.py` (sklearn DBSCAN, see above) — not `gmx cluster`. The PBC scripts
> are deliberately *not* named `*cluster*` to keep that distinction clear.

The multi-chain scripts (mirror their single-chain counterparts):

| Script | Mirrors | Role |
|--------|---------|------|
| `multichain_fix_PBC.sh` | `fix_PBC.sh` | PBC fix in 3 passes: `whole → cluster → mol+center+compact` |
| `multichain_fix_PBC_strip_align.sh` | `fix_PBC_strip_align.sh` | orchestrator (reuses `strip_and_align_trajectory.sh`) |
| `multichain_extract_protein.sh` | `extract_protein.sh` | RMSD reference via `-pbc cluster` |
| `multichain_chain_index.py` | — | chain detector + index builder (per-chain + per-chain-backbone groups) |
| `multichain_chain_rmsd.sh` | `calc_traj_rmsd.sh` | one chain's backbone RMSD, fit to itself |
| `multichain_chain_rmsf.sh` | `calc_traj_rmsf.sh` | one chain's per-residue RMSF |
| `multichain_interchain_dist.sh` | — | inter-chain minimum-image distance (`gmx mindist`) |

**Output mirrors the single-chain set** (`<p>_stripped_aligned.{xtc,gro}`,
`<p>_init.gro`, `<p>_{rmsd,rg,rmsf,dssp}.*` on the whole complex), **plus**:
`<p>_chain{A,B,…}_{rmsd,rmsf}.{xvg,png}` and `<p>_interchain_mindist.{xvg,png}`
(one curve per chain pair). Per-chain RMSD fits each chain to *itself*, so it reports
that chain's internal drift regardless of how the chains sit relative to each other.

`multichain_chain_index.py` uses the protein **molecule types** from the topology
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
XTC=OUTDIR/trajectories/md.xtc
P=OUTDIR/analysis/md

bash   fix_PBC_strip_align.sh "$TPR" "$XTC" "$P"
bash   calc_traj_rmsd.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rmsd.xvg" && python plot_xvg.py "${P}_rmsd.xvg" "${P}_rmsd.png"
bash   calc_traj_rg.sh   "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rg.xvg"   && python plot_xvg.py "${P}_rg.xvg"   "${P}_rg.png"
bash   calc_traj_rmsf.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_rmsf.xvg" && python plot_xvg.py "${P}_rmsf.xvg" "${P}_rmsf.png"
bash   calc_traj_dssp.sh "${P}_stripped_aligned.gro" "${P}_stripped_aligned.xtc" "${P}_dssp.dat" && python plot_dssp.py "${P}_dssp.dat" "${P}_dssp.png"
```

### T-REMD
Same as above, but the inputs are the lowest-T slot and the prefix is
`remd_rep000`, plus the REMD-only acceptance-rate QC:
```bash
TPR=OUTDIR/prod/rep000/remd.tpr
XTC=OUTDIR/trajectories/remd_rep000.xtc
P=OUTDIR/analysis/remd_rep000

python remd_acceptance.py OUTDIR              # QC: exchange acceptance rates
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
