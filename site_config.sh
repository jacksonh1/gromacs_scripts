#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# site_config.sh — cluster-level settings for the GROMACS REMD pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Edit this file once. It is sourced automatically by the engine scripts
# (REMD-gromacs.sbatch, etc.) via the GROMACS_SCRIPTS_DIR environment variable
# set in each submit script. You never need to touch it again after setup.
# ─────────────────────────────────────────────────────────────────────────────

# ── GROMACS installation (T-REMD / plain MD) ──────────────────────────────────
# Point this at the GROMACS bin/GMXRC used by REMD-gromacs.sbatch and MD-gromacs.sbatch.
GMXRC="${GMXRC:-$HOME/opt/gromacs/2024.3-plumed/bin/GMXRC}"

# ── REST2 (PLUMED Hamiltonian replica exchange) ───────────────────────────────
# REST2 needs a GROMACS build with WORKING PLUMED hrex. IMPORTANT: the 2024.3-plumed
# build above does NOT work for hrex — it runs but silently gives ZERO exchanges
# (hrex was never ported to the GROMACS-2024 PLUMED patch; see CLAUDE.md). REST2 uses
# a dedicated GROMACS 2023.5 build instead (scripts/installation/install_gromacs-2023.5-plumed.sh).
REST2_GMXRC="${REST2_GMXRC:-$HOME/opt/gromacs/2023.5-plumed/bin/GMXRC}"

# ── Custom force fields (GMXLIB) ──────────────────────────────────────────────
# Extra directory GROMACS searches for <name>.ff force-field dirs, IN ADDITION to
# each build's bundled share/gromacs/top (so amber99sb-ildn etc. keep working).
# This is how FF=charmm36m is found — the CHARMM36m GROMACS port lives here as
# charmm36m.ff (installed from the MacKerell force-switch port; -water tip3p then
# resolves to the CHARMM-modified TIP3P inside that dir). Covers both builds above.
# GMXRC does not set GMXLIB, so exporting it here survives sourcing GMXRC.
export GMXLIB="$HOME/opt/gromacs/ff${GMXLIB:+:$GMXLIB}"
# CUDA matching that build (2023.5 predates cuda 12.9; use 12.4 via deprecated-modules).
REST2_CUDA_MODULE="${REST2_CUDA_MODULE:-cuda/12.4.0}"

# Script that activates the PLUMED kernel (sets PLUMED_KERNEL, LD_LIBRARY_PATH). REST2 only.
PLUMED_SH="${PLUMED_SH:-$HOME/plumed.sh}"

# hcoll/ocoms compat libs. `gmx_mpi` AND the `plumed` CLI link libhcoll.so.1 / libocoms.so.0,
# which some GPU compute nodes lack under /opt/mellanox. The REST2 engine prepends this dir
# to LD_LIBRARY_PATH so both run everywhere (harmless no-op on nodes that already have them).
# Populate once:  mkdir -p ~/opt/hcoll_compat && cp -a /opt/mellanox/hcoll/lib/libhcoll.so.1* \
#                    /opt/mellanox/hcoll/lib/libocoms.so.0* ~/opt/hcoll_compat/   (from a node that HAS them)
HCOLL_COMPAT_DIR="${HCOLL_COMPAT_DIR:-$HOME/opt/hcoll_compat}"

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
