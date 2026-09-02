# GROMACS Plain MD Output Guide

Generated from `MD-gromacs.sbatch`. The top-level output directory is set by
`OUTDIR` in the submit script, e.g.:

```
outputs/output_MD/helix_fusion-100ns-MD-300K-NPT/
```

This is the single-temperature production pipeline (no replica exchange). It starts from
a designed/folded input pose — never an unfolded state — and is mainly used for
**bound-state ensemble sampling** of a complex in its bound pose (objective #4), and
optionally for single-structure **stability/rigidity** (#1) and **flexible-region**
identification (#2). Protocol:
**EM → heat (NVT, restrained) → density (NPT, restrained, iterative) →
[optional relax (NPT, unrestrained)] → production (NPT) → post-analysis.**

Position restraints (`-DPOSRES`) are held through EM/heat/density to preserve the
input structure. By **default (`RELAX_NS=0`)** the restraints are released at the
**start of production**, so the analyzed trajectory captures the protein relaxing away
from the designed pose — the signal for stability / flexibility / variant-comparison.
Set **`RELAX_NS>0`** to insert an unrestrained equilibration first, so production
starts pre-relaxed (e.g. for bound-state equilibrium sampling, where you'd rather drop
the relaxation transient up front than trim it later). Either way equilibration and
production share one ensemble (V-rescale thermostat + C-rescale barostat), so there is
no thermostat/barostat switch.

---

## Directory Tree

```
OUTDIR/
├── build/                  # Step 2 — System building                    [REAL]
├── em/                     # Step 3 — Energy minimization                 [REAL]
├── analysis/               # Step 10 — post-analysis outputs              [REAL]
├── solvated_snapshots/     # Step 10 — PDB snapshots incl. 5 Å solvent shell [REAL]
├── logs/                   # mdrun stdout logs                            [REAL]
├── heat/     ─▶ scratch    # Step 4 — Heat (NVT, restrained)              [SYMLINK]
├── density/  ─▶ scratch    # Step 5 — Density equilibration (NPT, restrained, iterative)
├── relax/    ─▶ scratch    # Step 6 — Relax (unrestrained NPT) — only if RELAX_NS>0
├── prod/     ─▶ scratch    # Step 7 — Production MD (NPT, unrestrained)   [SYMLINK]
│                           #   prod/md.xtc = the production trajectory (analyze directly)
├── parameters.txt          # Summary of all job parameters
└── <OUTBASE>_final.pdb     # Final structure exported from production
```

**Output model (folder-symlink).** The small dirs (`build/ em/ analysis/
solvated_snapshots/ logs/`,
`parameters.txt`, final PDB) are **real** in `OUTDIR`; the bulk stage dirs (`heat/
density/ [relax/] prod/`) are **folder symlinks into scratch** (`SCRATCH_DIR`), so
mdrun writes trajectories straight onto scratch — quota-safe, laptop-sync-friendly, and
analysis reads `prod/md.xtc` through the symlink with no path change. On success the
run copies the real dirs into `SCRATCH_DIR` so the scratch archive stands alone. Set
**`SYMLINK_BULK=0`** to keep the stage dirs **real in `OUTDIR`** with no scratch
offload; everything else here is identical.

---

## Step-by-step File Reference

### `build/` — System Building (Step 2)

Identical to the REMD pipeline. Converts the input PDB into a solvated, ionized
system. Canonical outputs: `<OUTBASE>_ions.gro` (final system) and
`<OUTBASE>.top` (topology), used by all downstream steps.

### `em/` — Energy Minimization (Step 3)

Steepest-descent minimization with position restraints on. `em.gro` is the
minimized structure, input to the heat step.

### `heat/` — Heat (NVT, Step 4)

Heats the system to `T_SIM` at constant volume, **with position restraints on**
the protein heavy atoms (`-DPOSRES`). Velocities are generated here
(`gen-vel = yes`).

| File | Description |
|------|-------------|
| `heat.mdp` | NVT parameters (V-rescale, `ref-t = T_SIM`, `gen-vel = yes`, `-DPOSRES`) |
| `heat.tpr` | Run input |
| `heat.log` | mdrun log |
| `heat.edr` | Energy file |
| `heat.gro` | Equilibrated structure — input to NPT |
| `heat.cpt` | Checkpoint — carries velocities forward into the density stage |

> `logs/mdrun_heat.log` contains the mdrun stdout.

### `density/` — Density Equilibration (NPT, Step 5)

Equilibrates density at `T_SIM` under constant pressure, **restraints still on**.
Continues from the heat checkpoint (`gen-vel = no`, `continuation = yes`), and every
segment thereafter resumes from the previous segment's checkpoint (`grompp -t`), so the
V-rescale and C-rescale internal state and RNG streams carry across each boundary.

Position restraints are anchored to the **minimized** structure (`em/em.gro`) for every
segment. A moving reference (the previous segment's output) would resist per-segment
displacement but never cumulative drift away from the input pose.

Runs iteratively, up to `DENSITY_MAX_SEG` segments, until the volume plateaus: a slope
test over the trailing `DENSITY_MIN_SEG` segments, not a last-two-segments comparison —
see `docs/PARAMETERS.md`.

| File | Description |
|------|-------------|
| `density_seg<N>.mdp` | Parameters for segment N (V-rescale + C-rescale, `-DPOSRES`) |
| `density_seg<N>.tpr` | Run input for segment N |
| `density_seg<N>.log` | mdrun log |
| `density_seg<N>.edr` | Energy file — volume extracted for convergence check |
| `density_seg<N>.gro` | Structure at end of segment N |
| `density_seg<N>.cpt` | Checkpoint |
| `volume_seg<N>.xvg` | Average volume from the edr (convergence check) |

> Trajectories (`.xtc`) go to scratch. The last converged `density_seg<N>.gro` +
> `.cpt` are the starting point for the unrestrained equilibration.

### `relax/` — Relax (unrestrained NPT, Step 6, optional)

**Only created when `RELAX_NS>0` (default `0` = skipped).** A short NPT run with
**restraints released** and the **production ensemble** (V-rescale + C-rescale),
continuing from the last NPT-density checkpoint. It lets the protein relax from its
restrained (input) pose and the barostat settle, so production starts pre-relaxed
instead of capturing the restraint-release transient — useful when you want production
to be an equilibrium ensemble (e.g. bound-state sampling). When skipped (the default),
production continues straight from the last restrained NPT-density segment, restraints
come off at production start, and the relaxation is recorded in the analyzed trajectory.

| File | Description |
|------|-------------|
| `relax.mdp` | Parameters (V-rescale + C-rescale, **no** `-DPOSRES`, `continuation = yes`) |
| `relax.tpr` | Run input |
| `relax.log` | mdrun log |
| `relax.edr` | Energy file |
| `relax.gro` | Equilibrated structure — input to production |
| `relax.cpt` | Checkpoint — carries velocities/box into production |

> `logs/mdrun_relax.log` contains mdrun stdout. The `.xtc` goes to scratch.

### `prod/` — Production MD (Step 7)

The production run: **NPT ensemble, no restraints**, V-rescale thermostat + C-rescale
barostat. Continues from the unrestrained `eq` checkpoint if that step ran, otherwise
straight from the last restrained NPT-density segment (in which case **restraints are
released here**, at production start). Coordinates, velocities, and box carry over
either way; the ensemble matches equilibration, so there is no thermostat/barostat switch.

| File | Description |
|------|-------------|
| `md.mdp` | Production parameters (V-rescale + C-rescale, no `-DPOSRES`) |
| `md.tpr` | Run input |
| `md.log` | mdrun log |
| `md.edr` | Energy file |
| `md.gro` | Final structure |
| `md.cpt` | Checkpoint (use to restart/extend) |

> `logs/mdrun_md.log` contains mdrun stdout. The production trajectory is
> `prod/md.xtc`, which lives on scratch via the `prod/` folder symlink (or in
> `OUTDIR` when `SYMLINK_BULK=0`).

### Where the trajectories live

There is **no `trajectories/` collection dir**. Under the folder-symlink model each
`.xtc` sits in its stage dir, which is a symlink into scratch (`SCRATCH_DIR`):
- `prod/md.xtc` — production trajectory (the analysis target)
- `heat/heat.xtc` — heat-stage (NVT) trajectory
- `relax/relax.xtc` — relax-stage (unrestrained NPT) trajectory (only if `RELAX_NS>0`)
- `density/<OUTBASE>_density_seg<N>.xtc` — NPT density equilibration trajectories

> These resolve to `SCRATCH_DIR` under `SCRATCH_ROOT`, which also holds a self-contained
> copy of `analysis/ logs/ em/ build/` + `parameters.txt`. Copy anything you need
> long-term. With `SYMLINK_BULK=0` the stage dirs are real in `OUTDIR` (no scratch copy).

### `analysis/` — Post-analysis (Step 10)

Produced by the shared `scripts/analysis/` tools (see `scripts/analysis/README.md`):

| File | Description |
|------|-------------|
| `md_stripped_aligned.xtc` | Protein-only, PBC-fixed, backbone-aligned trajectory |
| `md_stripped_aligned.gro` | Protein-only first-frame reference (used for Rg/RMSF/DSSP) |
| `md_init.gro` | Protein-only **initial structure** (from minimized `em.gro`, made whole with `-pbc whole`); the RMSD reference |
| `md_rmsd.xvg` / `.png` | Backbone RMSD vs time — drift from the initial structure (`md_init.gro`). Plot in Å; `.xvg` in nm |
| `md_rg.xvg` / `.png` | Radius of gyration vs time |
| `md_rmsf.xvg` / `.png` | Per-residue backbone RMSF |
| `md_dssp.dat` / `.png` | DSSP secondary structure (per-frame data + residue×frame map) |
| `clustering/md_cluster_assignments.csv` | Per-frame conformational-cluster assignment (`frame,time_ps,cluster`) |
| `clustering/md_cluster_rep_c00.pdb`, `_c01.pdb`, … | Representative structure per cluster (ranked by population) |
| `clustering/md_cluster_populations.png` / `_timeseries.png` | Cluster populations + cluster-vs-time |
| `clustering/md_cluster_summary.txt` | Clustering method/cutoff/min_samples + per-cluster population table |

**Multi-chain complexes** (auto-detected from the topology): the analysis uses the
`multichain_*` pipeline (keeps the chains in one periodic image via `-pbc cluster`),
produces the **same files as above** on the whole complex, **plus**:

| File | Description |
|------|-------------|
| `md_chain{A,B,…}_rmsd.xvg` / `.png` | per-chain backbone RMSD, each chain fit to itself |
| `md_chain{A,B,…}_rmsf.xvg` / `.png` | per-chain per-residue RMSF |
| `md_interchain_mindist.xvg` / `.png` | minimum-image distance between chains (binding check) |
| `md_chains.ndx` | per-chain index groups used for the above |

See `scripts/analysis/README.md` → "Multi-chain complexes".

### `solvated_snapshots/` — PDBs that keep the solvent (Step 10, MD only)

Everything in `analysis/` is protein-only. These are the exception: a handful of
PBC-corrected, mutually-aligned PDBs holding the protein **plus every whole water/ion
within 5 Å of it** — for inspecting interface / bound-state waters.

| File | Description |
|------|-------------|
| `md_solvshell_<NNNNNN>ps.pdb` | One snapshot, protein + 5 Å solvent shell, backbone-fitted to the t=0 snapshot. Zero-padded ps ⇒ lexical order is chronological |

Five snapshots by default (first, three middle, last). Load them all at once — they
superimpose. Both the shell thickness and the count are overridable:
`SHELL_NM=0.8 N_SNAPSHOTS=9 bash scripts/analysis/run_analysis.sh OUTDIR`.

The atom count differs between snapshots **by design** — the shell is recomputed per
frame, so this is not a trajectory and there is no matching topology file. See
`scripts/analysis/README.md` → "Solvated snapshots" for why.

### `parameters.txt`

Plain-text record of all simulation parameters (force field, temperature,
timestep, production length, NPT convergence segments, etc.) and the scratch path.

It ends with a **Provenance** block recording what actually ran, so a finished job
is self-describing without the submit script: the engine path, the resolved
`site_config.sh`, the GROMACS version, `GMXRC`, `GMXLIB` (which is how
a `FF=charmm*` name resolves), the input structure, and `$SLURM_SUBMIT_DIR`.

---

## Key Files for Analysis

| Goal | File(s) |
|------|---------|
| Check EM converged | `em/em.log` — look for `Fmax <` line |
| Check density equilibration | `density/volume_seg*.xvg` |
| Analyze the ensemble | `analysis/md_stripped_aligned.xtc` (protein-only, aligned) + `analysis/md_stripped_aligned.gro` |
| RMSD / Rg / RMSF / DSSP | `analysis/md_{rmsd,rg,rmsf,dssp}.*` |
| Conformational states | `analysis/clustering/md_cluster_summary.txt` + `md_cluster_rep_c*.pdb` (representative structures) |
| Re-run post-processing | `bash scripts/analysis/run_analysis.sh OUTDIR` (regenerates the whole `analysis/` dir; no resubmission) |
| Inspect interface waters | `solvated_snapshots/md_solvshell_*.pdb` — open all 5 at once (already aligned) |
| Inspect with full solvent | `bash scripts/analysis/fix_PBC.sh prod/md.tpr prod/md.xtc analysis/md_pbc.xtc` |
| Final structure | `<OUTBASE>_final.pdb` — protein-centred, compact cell (`-pbc mol -center -ur compact`) |
