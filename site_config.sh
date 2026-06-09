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

# ── Python / conda environment (analysis + plotting) ──────────────────────────
# The post-analysis tools (matplotlib, mdanalysis, numpy, …) run in a conda env
# created from scripts/installation/environment.yml. Create it once with
# scripts/installation/install_python_env.sh.
#
# CONDA_MODULE — module that provides conda/mamba (miniforge on this cluster).
# GROMD_ENV    — name of the conda env (matches `name:` in environment.yml).
CONDA_MODULE="${CONDA_MODULE:-miniforge/25.11.0-0}"
GROMD_ENV="${GROMD_ENV:-groMD_env}"

# Activate the analysis/plotting conda env. Call this immediately before the
# Python analysis tools — NOT around the mdrun steps, so conda's libraries can't
# interfere with the GROMACS MPI/CUDA runtime. Preserves the caller's `set -u`
# state (conda's activation scripts reference unset variables).
activate_python_env() {
  module load "$CONDA_MODULE"
  local had_u=0; [[ $- == *u* ]] && had_u=1
  set +u
  eval "$(conda shell.bash hook)"
  conda activate "$GROMD_ENV"
  (( had_u )) && set -u
  return 0
}
