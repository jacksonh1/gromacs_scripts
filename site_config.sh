#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# site_config.sh — cluster-level settings for the GROMACS REMD pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Edit this file once. It is sourced automatically by the engine scripts
# (REMD-gromacs.sbatch, etc.) via the GROMACS_SCRIPTS_DIR environment variable
# set in each submit script. You never need to touch it again after setup.
# ─────────────────────────────────────────────────────────────────────────────

# ── GROMACS installation ──────────────────────────────────────────────────────
# Point this at your GROMACS bin/GMXRC. A PLUMED-patched build works for all
# standard T-REMD tasks and additionally enables the REST2 pipeline. A plain
# build works for T-REMD only (REST2 will fail with a clear error).
GMXRC="${GMXRC:-$HOME/opt/gromacs/2024.3-plumed/bin/GMXRC}"

# Script that activates the PLUMED kernel (sets PLUMED_KERNEL, LD_LIBRARY_PATH).
# Only needed for REST2; ignored by the T-REMD pipeline.
PLUMED_SH="${PLUMED_SH:-$HOME/plumed.sh}"

# ── Scratch storage ───────────────────────────────────────────────────────────
# Root directory for large trajectory files. Each job creates a timestamped
# subdirectory here. Requires ~100 GB free per job.
SCRATCH_ROOT="${SCRATCH_ROOT:-/orcd/data/keating/001/${USER}/MD}"

# ── Cluster modules ───────────────────────────────────────────────────────────
CUDA_MODULE="${CUDA_MODULE:-cuda/12.9.1}"
OPENMPI_MODULE="${OPENMPI_MODULE:-openmpi/5.0.8}"
