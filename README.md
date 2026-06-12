# GROMACS Structure-Characterization Pipelines - UNSTABLE

# warning - This pipeline is under active development


GROMACS 2024.3 pipelines for characterizing **designed protein structures** on single-node GPU clusters (SLURM). Built for the Keating lab at MIT; configurable for any cluster via `site_config.sh`.

Two engines are provided:
- **T-REMD** (`submit_REMD.sh`) — temperature replica-exchange enhanced sampling
- **Plain production MD** (`submit_MD.sh`) — single-temperature NPT production

## Purpose

Every run starts from a folded/designed input structure — **never** from an unfolded or extended state. The tools serve four objectives:

1. **Stability/rigidity** characterization of a given starting (designed) structure
2. Identifying **flexible regions** (per-residue)
3. **Variant comparison** — run the same protocol on several design variants and see which best retains its designed conformation
4. **Bound-state sampling** — simulate a complex in its bound pose and sample the bound ensemble

Plain MD is mainly for **#4** (and optionally #1, #2); T-REMD is primarily for **#1–3**. The input pose is the reference the analysis is measured against (RMSD = drift from the design, RMSF = local flexibility).

---

## Prerequisites

- GROMACS 2024.3 compiled with MPI + CUDA (see `scripts/installation/`)
- A conda env for post-analysis (matplotlib, mdanalysis, numpy, …), created from `scripts/installation/environment.yml` — see One-Time Setup. The sbatch engines use the system `python3` (standard library only) for inline temperature-ladder and convergence calculations, and activate the conda env only for the post-analysis/plotting steps.
- SLURM with GPU access

A PLUMED 2.9.4-patched build is recommended — it runs all T-REMD tasks identically to a plain build and additionally enables the REST2 pipeline (`dev/`). A plain build works for T-REMD only.

---

## One-Time Setup

1. **Edit `site_config.sh`** in the repo root. Set:
   - `GMXRC` — path to your `bin/GMXRC` from the GROMACS installation
   - `SCRATCH_ROOT` — fast scratch storage root for trajectory files (~100 GB per job)
   - `CUDA_MODULE`, `OPENMPI_MODULE` — your cluster's module names
   - `CONDA_MODULE`, `GROMD_ENV` — module providing conda and the analysis env name

2. **Create the analysis conda env** (once per cluster, on a login node):
   ```bash
   bash scripts/installation/install_python_env.sh
   ```
   This builds `groMD_env` from `scripts/installation/environment.yml`. The
   pipeline activates it automatically for the post-analysis steps.

3. That's it. The pipeline scripts source `site_config.sh` automatically.

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

Parameters are validated at job start: a missing `PDB_IN`, an out-of-range or
non-numeric value, or (in a config file) a misspelled key fails the job immediately with
a clear `[ERROR]` message instead of silently falling back to a default.

---

## Pipeline Overview

The engine script (`scripts/simulation/REMD-gromacs.sbatch`) runs 13 steps:

| Step | Description |
|------|-------------|
| 0 | Load environment (modules, GROMACS) |
| 1 | Set parameters, create scratch and output directories |
| 2 | Compute geometric temperature ladder |
| 3 | Build system: pdb2gmx → editconf → solvate → genion |
| 4 | Energy minimization (steepest descent) |
| 5 | NPT density equilibration (iterative, convergence-checked) |
| 6 | Prepare per-replica equilibration inputs (NVT, or NPT if `ENSEMBLE=NPT`) |
| 7 | Run per-replica equilibration (all replicas in parallel via MPI) |
| 8 | Prepare REMD production inputs |
| 9 | Run T-REMD production |
| 10 | Finalize outputs, create trajectory symlinks |
| 11 | Write parameters log |
| 12 | Post-analysis: acceptance rates + PBC/strip/align + RMSD/Rg/RMSF/DSSP + clustering (rep000) |

See `scripts/simulation/REMD-output-guide.md` for a full description of all output files.

### Production ensemble (`ENSEMBLE=NVT|NPT`)

T-REMD production runs **NVT** by default (constant volume). Set `ENSEMBLE=NPT` in the
submit script to pressure-couple production with a **C-rescale barostat** (at `REF_P` bar,
default 1.0) — the correct ensemble for density-sensitive observables and bound-state
sampling. Under `NPT` the per-replica `equil/` stage is also pressure-coupled, and GROMACS
automatically adds the *PV* term to the replica-exchange Metropolis criterion (`-replex` is
unchanged). The constant-temperature interpretation is unaffected either way:
`prod/rep000/remd.xtc` is still the lowest-temperature ensemble. Plain MD production is
always NPT.

---

## Output & Analysis

Large trajectory files (`.xtc`) are written to `SCRATCH_ROOT` and symlinked into `OUTDIR/trajectories/`. Copy them before scratch is purged.

**Key point:** In GROMACS T-REMD, each replica runs at a **fixed temperature** and **coordinates** are exchanged between replicas. This means:

- `prod/rep000/remd.xtc` is the 300 K constant-temperature trajectory — use it directly for analysis.
- No demux step is needed.

Post-analysis runs automatically at the end of each job. To re-run manually:
```bash
python scripts/analysis/remd_acceptance.py    OUTDIR        # acceptance rates (target 20–30%)
bash   scripts/analysis/fix_PBC_strip_align.sh OUTDIR 000   # PBC fix + strip + align
```
See `scripts/analysis/README.md` for the full script reference.

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
├── scripts/
│   ├── simulation/             # Simulation engines (do not edit)
│   │   ├── REMD-gromacs.sbatch    # T-REMD engine
│   │   ├── MD-gromacs.sbatch      # Plain-MD engine
│   │   ├── config_example.sh      # Job config template (copy and edit)
│   │   └── REMD-output-guide.md   # Full output file reference
│   ├── analysis/               # Post-processing tools (see scripts/analysis/README.md)
│   └── installation/           # GROMACS + PLUMED build scripts
├── dev/                        # REST2 pipeline (in development)
└── example/
    ├── input_pdbs/             # Example protein structures
    └── submit_jobs/
        └── submit_REMD.sh     # Submission wrapper (reads site_config.sh)
```
