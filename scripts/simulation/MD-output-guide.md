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
├── build/                  # Step 2 — System building
├── em/                     # Step 3 — Energy minimization
├── heat/                   # Step 4 — Heat (NVT, restrained)
├── density/                # Step 5 — Density equilibration (NPT, restrained, iterative)
├── relax/                  # Step 6 — Relax (unrestrained NPT) — only if RELAX_NS>0
├── prod/                   # Step 7 — Production MD (NPT, unrestrained)
├── logs/                   # mdrun stdout logs
├── trajectories/           # Symlinks to scratch trajectories
├── analysis/               # Step 10 — post-analysis outputs
├── parameters.txt          # Summary of all job parameters
└── <OUTBASE>_final.pdb     # Final structure exported from production
```

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
Continues from the heat checkpoint (`gen-vel = no`, `continuation = yes`). Runs
iteratively (up to `DENSITY_MAX_SEG` segments) until volume converges
(`< DENSITY_TOL_REL` relative change over `DENSITY_MIN_SEG`+ segments).

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

> `logs/mdrun_md.log` contains mdrun stdout. The trajectory (`md.xtc`) is written
> to scratch and symlinked into `trajectories/`.

### `trajectories/`

Symlinks pointing to the actual trajectory files on scratch (`SCRATCH_DIR`):
- `md.xtc` — production trajectory
- `heat.xtc` — heat-stage (NVT) trajectory
- `relax.xtc` — relax-stage (unrestrained NPT) trajectory (only if `RELAX_NS>0`)
- `<OUTBASE>_density_seg<N>.xtc` — NPT density equilibration trajectories

> These live under `SCRATCH_ROOT`. Copy them before scratch is purged.

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

### `parameters.txt`

Plain-text record of all simulation parameters (force field, temperature,
timestep, production length, NPT convergence segments, etc.) and the scratch path.

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
| Inspect with solvent | `bash scripts/analysis/fix_PBC.sh prod/md.tpr trajectories/md.xtc analysis/md_pbc.xtc` |
