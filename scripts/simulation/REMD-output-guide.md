# GROMACS T-REMD Output Guide

Temperature replica-exchange enhanced sampling, starting from a designed/folded input
pose — never an unfolded state. Primarily used for **stability/rigidity** characterization
(objective #1), **flexible-region** identification (#2), and **variant comparison** —
running the same protocol on several design variants to see which best retains its
designed conformation (#3).

Generated from `REMD-gromacs.sbatch`. The top-level output directory is set by `OUTDIR` in the submit script, e.g.:

```
outputs/output_T-REMD-gromacs/DB2_unbound-5ns-REMD-300-450K-48reps-NVT/
```

---

## Directory Tree

```
OUTDIR/
├── build/                        # Step 3 — System building
├── em/                           # Step 4 — Energy minimization
├── density/                          # Step 5 — NPT density equilibration
├── equil/                          # Steps 6–7 — Per-replica equilibration (NVT, or NPT if ENSEMBLE=NPT)
│   ├── rep000/
│   ├── rep001/
│   └── ... (one per replica)
├── prod/                         # Steps 8–9 — REMD production
│   ├── rep000/
│   ├── rep001/
│   └── ... (one per replica)
├── logs/                         # mdrun stdout logs
├── trajectories/                 # Symlinks to scratch trajectories
├── parameters.txt                # Summary of all job parameters
└── <OUTBASE>_final_rep000.pdb    # Final structure exported from replica 000
```

---

## Step-by-step File Reference

### `build/` — System Building (Step 3)

Converts the input PDB into a solvated, ionized GROMACS system.

| File | Produced by | Description |
|------|-------------|-------------|
| `pdb2gmx.log` | `pdb2gmx` | Log for topology generation |
| `<OUTBASE>_processed.gro` | `pdb2gmx` | Structure with H atoms added, ready for GROMACS |
| `<OUTBASE>.top` | `pdb2gmx` | Master topology file (includes ff parameters, water, ions) |
| `<OUTBASE>_posre.itp` | `pdb2gmx` | Position restraint definitions for protein heavy atoms |
| `editconf.log` | `editconf` | Log for box setup |
| `<OUTBASE>_box.gro` | `editconf` | Structure with simulation box defined (dodecahedron by default) |
| `solvate.log` | `solvate` | Log for water addition |
| `<OUTBASE>_solv.gro` | `solvate` | Solvated structure |
| `ions.mdp` | script | Minimal mdp used only to generate a .tpr for `genion` |
| `grompp_ions.log` | `grompp` | Log for ion .tpr generation |
| `ions.tpr` | `grompp` | Run input for `genion` (not a real simulation) |
| `genion.log` | `genion` | Log for ion placement |
| `<OUTBASE>_ions.gro` | `genion` | **Final system**: solvated + neutralized + 0.15 M NaCl |

> `<OUTBASE>_ions.gro` and `<OUTBASE>.top` are the canonical system files used by all downstream steps.

---

### `em/` — Energy Minimization (Step 4)

Relaxes clashes and strained geometry from system building. No dynamics — atoms are moved purely to reduce potential energy.

| File | Description |
|------|-------------|
| `em.mdp` | EM parameters (steepest descent, 3000 steps, position restraints on) |
| `grompp_em.log` | `grompp` preprocessing log |
| `em.tpr` | Binary run input for EM |
| `em.log` | GROMACS mdrun log (force convergence, potential energy) |
| `em.edr` | Energy file (readable with `gmx energy`) |
| `em.gro` | **Minimized structure** — input to density equilibration |
| `em.trr` | Full-precision trajectory (usually not needed) |
| `mdrun_em.log` | mdrun stdout/stderr |

> Check `em.log` to verify the system reached the force tolerance (`Fmax < 1000 kJ/mol/nm`).

---

### `density/` — NPT Density Equilibration (Step 5)

Equilibrates system density at the lowest replica temperature (`T_MIN`) under constant pressure. Runs iteratively (up to `DENSITY_MAX_SEG` segments) until volume converges.

| File | Description |
|------|-------------|
| `density_seg<N>.mdp` | Parameters for segment N (V-rescale thermostat, C-rescale barostat) |
| `grompp_density<N>.log` | `grompp` log for segment N |
| `density_seg<N>.tpr` | Run input for segment N |
| `density_seg<N>.log` | GROMACS mdrun log |
| `density_seg<N>.edr` | Energy file — used to extract volume for convergence check |
| `density_seg<N>.gro` | Structure at end of segment N |
| `density_seg<N>.cpt` | Checkpoint file |
| `volume_seg<N>.xvg` | Average volume extracted from the edr (convergence check) |
| `mdrun_density<N>.log` | mdrun stdout/stderr |

> Trajectories (`.xtc`) go to scratch (`SCRATCH_DIR`), not here.
> The last converged `density_seg<N>.gro` is used as the starting structure for all replicas.

---

### `equil/rep<NNN>/` — Per-Replica Equilibration (Steps 6–7)

Each replica equilibrates at its own target temperature. Velocities are freshly generated for each replica. Constant volume by default (`ENSEMBLE=NVT`); with `ENSEMBLE=NPT` this stage also runs the C-rescale barostat so each replica relaxes its box at its own temperature before production.

| File | Description |
|------|-------------|
| `equil.mdp` | Equilibration parameters with `ref-t` set to this replica's temperature (NVT, or NPT C-rescale when `ENSEMBLE=NPT`) |
| `grompp.log` | `grompp` log |
| `equil.tpr` | Run input |
| `equil.log` | GROMACS mdrun log |
| `equil.edr` | Energy file |
| `equil.gro` | **Equilibrated structure** at this replica's T — input to production |
| `equil.cpt` | Checkpoint — carries velocities forward into production |

> `logs/mdrun_equil.log` contains the combined stdout from the MPI mdrun across all replicas.

---

### `prod/rep<NNN>/` — REMD Production (Steps 8–9)

The actual T-REMD simulation. Replicas run simultaneously and attempt coordinate exchanges every `REPLEX_PS` ps.

**Ensemble (`ENSEMBLE`, default `NVT`):** production runs at constant volume (`NVT`) or constant pressure (`NPT`, C-rescale barostat at `REF_P` bar). The constant-temperature interpretation is unchanged either way — `rep000/remd.xtc` is the `T_MIN` ensemble. Under `NPT`, GROMACS automatically adds the *PV* term to the replica-exchange criterion (no `-replex` change), and the per-replica `equil/` stage also runs NPT so each replica's box is pre-relaxed.

| File | Description |
|------|-------------|
| `remd.mdp` | Production parameters (`continuation = yes`, no velocity generation) |
| `grompp.log` | `grompp` log |
| `remd.tpr` | Run input |
| `remd.log` | **Key log**: contains exchange attempt statistics and acceptance rates |
| `remd.edr` | Energy file |
| `remd.gro` | Final structure at end of production |
| `remd.cpt` | Checkpoint (use to restart/extend the run) |

> `logs/mdrun_remd.log` contains combined mdrun stdout.
> Trajectories (`.xtc`) are written to scratch and symlinked into `trajectories/`.

---

### `logs/`

| File | Description |
|------|-------------|
| `mdrun_equil.log` | Combined stdout from all-replica equilibration mdrun |
| `mdrun_remd.log` | Combined stdout from REMD production mdrun |

---

### `trajectories/`

Symlinks pointing to the actual trajectory files on scratch (`SCRATCH_DIR`). Named:
- `equil_rep<NNN>.xtc` — per-replica equilibration trajectory
- `remd_rep<NNN>.xtc` — production trajectory per replica
- `<OUTBASE>_density_seg<N>.xtc` — density equilibration trajectories

> These files live in `/orcd/data/keating/001/<user>/MD/<jobid>_<timestamp>/`. Copy them before the scratch is purged.

---

### `parameters.txt`

A plain-text record of all simulation parameters (force field, temperatures, timestep, exchange interval, etc.) and the scratch directory path. Written at the end of the job.

---

## Key Files for Analysis

| Goal | File(s) |
|------|---------|
| Check EM converged | `em/em.log` — look for `Fmax <` line |
| Check density equilibration | `density/volume_seg*.xvg` |
| Check REMD exchange rates | `python scripts/analysis/remd_acceptance.py OUTDIR` — outputs table + CSV |
| Analyze lowest-T ensemble | `analysis/remd_rep000_stripped_aligned.xtc` (protein-only, aligned) + `analysis/remd_rep000_stripped_aligned.gro` |
| RMSD / Rg / RMSF / DSSP | `analysis/remd_rep000_{rmsd,rg,rmsf,dssp}.*` |
| Conformational states | `analysis/clustering/remd_rep000_cluster_summary.txt` + `remd_rep000_cluster_rep_c*.pdb` (representative structures); `_cluster_{populations,timeseries}.png`, `_cluster_assignments.csv` |
| Re-run post-processing | `bash scripts/analysis/run_analysis.sh OUTDIR 000` (regenerates the whole `analysis/` dir; no resubmission) |
| Inspect with solvent | `bash scripts/analysis/fix_PBC.sh prod/rep000/remd.tpr trajectories/remd_rep000.xtc analysis/remd_rep000_pbc.xtc` — keeps full-system trajectory |

> **Note on demux:** `gmx demux` follows a single *configuration* as it walks through temperature space. It is NOT needed to obtain the constant-T ensemble. `rep000/remd.xtc` already is the 300 K ensemble — each replica runs at a fixed temperature and coordinates are exchanged between replicas, so the slot trajectory is the correct thermodynamic ensemble at that temperature.

> **Multi-chain complexes:** if the system has >1 protein chain, the analysis auto-uses the `multichain_*` pipeline (keeps the chains in one periodic image via `-pbc cluster`, REMD-safe) and adds per-chain RMSD/RMSF + an inter-chain minimum-distance curve (`remd_rep000_chain{A,B,…}_{rmsd,rmsf}.*`, `remd_rep000_interchain_mindist.*`). Core files are unchanged. See `scripts/analysis/README.md` → "Multi-chain complexes".
