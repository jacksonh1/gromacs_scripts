# GROMACS REMD Pipeline - UNSTABLE

# warning - This pipeline is under active development


Temperature Replica Exchange MD (T-REMD) pipeline for GROMACS 2024.3, designed for single-node GPU clusters (SLURM). Built for the Keating lab at MIT; configurable for any cluster via `site_config.sh`.

---

## Prerequisites

- GROMACS 2024.3 compiled with MPI + CUDA (see `installation_scripts/`)
- Python 3 (standard library only; used for inline temperature ladder and convergence calculations)
- SLURM with GPU access

A PLUMED 2.9.4-patched build is recommended — it runs all T-REMD tasks identically to a plain build and additionally enables the REST2 pipeline (`dev/`). A plain build works for T-REMD only.

---

## One-Time Setup

1. **Edit `site_config.sh`** in the repo root. Set:
   - `GMXRC` — path to your `bin/GMXRC` from the GROMACS installation
   - `SCRATCH_ROOT` — fast scratch storage root for trajectory files (~100 GB per job)
   - `SBATCH_PARTITION`, `SBATCH_GPU_TYPE`, `SBATCH_GPUS_PER_NODE` — your cluster's SLURM settings

2. That's it. The pipeline scripts source `site_config.sh` automatically.

---

## Running a Job

1. Copy the example submit script and edit your parameters directly in it:
   ```bash
   cp example/submit_jobs/submit_REMD.sh my_job.sh
   # edit my_job.sh: set PDB_IN, REPLICAS, T_MAX, TOTAL_NS, etc.
   ```

2. Run it:
   ```bash
   bash my_job.sh
   ```

All job parameters live in the submit script — no separate config file needed.

---

## Pipeline Overview

The engine script (`gromacs_scripts/REMD-gromacs.sbatch`) runs 12 steps:

| Step | Description |
|------|-------------|
| 0 | Load environment (modules, GROMACS) |
| 1 | Set parameters, create scratch and output directories |
| 2 | Compute geometric temperature ladder |
| 3 | Build system: pdb2gmx → editconf → solvate → genion |
| 4 | Energy minimization (steepest descent) |
| 5 | NPT density equilibration (iterative, convergence-checked) |
| 6 | Prepare per-replica NVT equilibration inputs |
| 7 | Run NVT equilibration (all replicas in parallel via MPI) |
| 8 | Prepare REMD production inputs |
| 9 | Run T-REMD production |
| 10 | Finalize outputs, create trajectory symlinks |
| 11 | Write parameters log |

See `gromacs_scripts/REMD-output-guide.md` for a full description of all output files.

---

## Output & Analysis

Large trajectory files (`.xtc`) are written to `SCRATCH_ROOT` and symlinked into `OUTDIR/trajectories/`. Copy them before scratch is purged.

**Key point:** In GROMACS T-REMD, each replica runs at a **fixed temperature** and **coordinates** are exchanged between replicas. This means:

- `prod/rep000/remd.xtc` is the 300 K constant-temperature trajectory — use it directly for analysis.
- No demux step is needed.

Check exchange acceptance rates:
```bash
grep 'Repl  average probabilities' OUTDIR/prod/rep000/remd.log
```
Target: ~20–40% acceptance between adjacent replicas.

---

## Configuration Reference

| File | Who edits it | What it controls |
|------|-------------|-----------------|
| `site_config.sh` | Once per user/cluster | GROMACS path, scratch root, module names |
| `my_job.sh` (copy of `submit_REMD.sh`) | Per job | PDB input, replicas, temperature range, simulation length |
| `#SBATCH` headers in engine script | Only if changing resource defaults | Partition, GPU type, memory, wall time |

---

## Directory Structure

```
gromacs_REMD/
├── site_config.sh              # Cluster-level settings (edit once)
├── gromacs_scripts/
│   ├── REMD-gromacs.sbatch    # T-REMD engine (do not edit)
│   ├── config_example.sh      # Job config template (copy and edit)
│   └── REMD-output-guide.md   # Full output file reference
├── analysis_scripts/           # Post-processing tools (to be added)
├── dev/                        # REST2 pipeline (in development)
├── installation_scripts/       # GROMACS + PLUMED build scripts
└── example/
    ├── input_pdbs/             # Example protein structures
    └── submit_jobs/
        └── submit_REMD.sh     # Submission wrapper (reads site_config.sh)
```
