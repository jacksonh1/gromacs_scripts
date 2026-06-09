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
├── npt/                          # Step 5 — NPT density equilibration
├── nvt/                          # Steps 6–7 — Per-replica NVT equilibration
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
| `em.gro` | **Minimized structure** — input to NPT equilibration |
| `em.trr` | Full-precision trajectory (usually not needed) |
| `mdrun_em.log` | mdrun stdout/stderr |

> Check `em.log` to verify the system reached the force tolerance (`Fmax < 1000 kJ/mol/nm`).

---

### `npt/` — NPT Density Equilibration (Step 5)

Equilibrates system density at the lowest replica temperature (`T_MIN`) under constant pressure. Runs iteratively (up to `NPT_MAX_SEG` segments) until volume converges.

| File | Description |
|------|-------------|
| `npt_seg<N>.mdp` | Parameters for segment N (V-rescale thermostat, C-rescale barostat) |
| `grompp_npt<N>.log` | `grompp` log for segment N |
| `npt_seg<N>.tpr` | Run input for segment N |
| `npt_seg<N>.log` | GROMACS mdrun log |
| `npt_seg<N>.edr` | Energy file — used to extract volume for convergence check |
| `npt_seg<N>.gro` | Structure at end of segment N |
| `npt_seg<N>.cpt` | Checkpoint file |
| `volume_seg<N>.xvg` | Average volume extracted from the edr (convergence check) |
| `mdrun_npt<N>.log` | mdrun stdout/stderr |

> Trajectories (`.xtc`) go to scratch (`SCRATCH_DIR`), not here.
> The last converged `npt_seg<N>.gro` is used as the starting structure for all replicas.

---

### `nvt/rep<NNN>/` — Per-Replica NVT Equilibration (Steps 6–7)

Each replica equilibrates at its own target temperature with no pressure coupling. Velocities are freshly generated for each replica.

| File | Description |
|------|-------------|
| `nvt.mdp` | NVT parameters with `ref-t` set to this replica's temperature |
| `grompp.log` | `grompp` log |
| `nvt.tpr` | Run input |
| `nvt.log` | GROMACS mdrun log |
| `nvt.edr` | Energy file |
| `nvt.gro` | **Equilibrated structure** at this replica's T — input to production |
| `nvt.cpt` | Checkpoint — carries velocities forward into production |

> `logs/mdrun_nvt_equil.log` contains the combined stdout from the MPI mdrun across all replicas.

---

### `prod/rep<NNN>/` — REMD Production (Steps 8–9)

The actual T-REMD simulation. Replicas run simultaneously and attempt coordinate exchanges every `REPLEX_PS` ps.

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
| `mdrun_nvt_equil.log` | Combined stdout from all-replica NVT mdrun |
| `mdrun_remd.log` | Combined stdout from REMD production mdrun |

---

### `trajectories/`

Symlinks pointing to the actual trajectory files on scratch (`SCRATCH_DIR`). Named:
- `nvt_rep<NNN>.xtc` — NVT equilibration trajectory per replica
- `remd_rep<NNN>.xtc` — production trajectory per replica
- `<OUTBASE>_npt_seg<N>.xtc` — NPT equilibration trajectories

> These files live in `/orcd/data/keating/001/<user>/MD/<jobid>_<timestamp>/`. Copy them before the scratch is purged.

---

### `parameters.txt`

A plain-text record of all simulation parameters (force field, temperatures, timestep, exchange interval, etc.) and the scratch directory path. Written at the end of the job.

---

## Key Files for Analysis

| Goal | File(s) |
|------|---------|
| Check EM converged | `em/em.log` — look for `Fmax <` line |
| Check density equilibration | `npt/volume_seg*.xvg` |
| Check REMD exchange rates | `python scripts/analysis/remd_acceptance.py OUTDIR` — outputs table + CSV |
| Analyze lowest-T ensemble | `analysis/remd_rep000_stripped_aligned.xtc` (protein-only, aligned) + `analysis/remd_rep000_stripped_aligned.gro` |
| RMSD / Rg / RMSF / DSSP | `analysis/remd_rep000_{rmsd,rg,rmsf,dssp}.*` |
| Conformational states | `analysis/clustering/remd_rep000_cluster_summary.txt` + `remd_rep000_cluster_rep_c*.pdb` (representative structures); `_cluster_{populations,timeseries}.png`, `_cluster_assignments.csv` |
| Re-run post-processing | `bash scripts/analysis/run_analysis.sh OUTDIR 000` (regenerates the whole `analysis/` dir; no resubmission) |
| Inspect with solvent | `bash scripts/analysis/fix_PBC.sh prod/rep000/remd.tpr trajectories/remd_rep000.xtc analysis/remd_rep000_pbc.xtc` — keeps full-system trajectory |

> **Note on demux:** `gmx demux` follows a single *configuration* as it walks through temperature space. It is NOT needed to obtain the constant-T ensemble. `rep000/remd.xtc` already is the 300 K ensemble — each replica runs at a fixed temperature and coordinates are exchanged between replicas, so the slot trajectory is the correct thermodynamic ensemble at that temperature.

> **Multi-chain complexes:** if the system has >1 protein chain, the analysis auto-uses the `multichain_*` pipeline (keeps the chains in one periodic image via `-pbc cluster`, REMD-safe) and adds per-chain RMSD/RMSF + an inter-chain minimum-distance curve (`remd_rep000_chain{A,B,…}_{rmsd,rmsf}.*`, `remd_rep000_interchain_mindist.*`). Core files are unchanged. See `scripts/analysis/README.md` → "Multi-chain complexes".
