#!/usr/bin/env bash
# Submit a T-REMD job.
# Usage: bash submit_REMD.sh path/to/job_config.sh
#
# job_config.sh sets simulation parameters (PDB_IN, REPLICAS, T_MAX, etc.).
# See scripts/simulation/config_example.sh for a template.
#
# Partition, GPU type, and other SLURM defaults are set by the #SBATCH headers
# in scripts/simulation/REMD-gromacs.sbatch. Override per-job with sbatch flags:
#   bash submit_REMD.sh my_config.sh -p other_partition --gres=gpu:a100:2

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

JOB_CONFIG="${1:-}"
shift || true   # consume $1 so remaining args (extra sbatch flags) are in "$@"

[[ -n "$JOB_CONFIG" && ! -f "$JOB_CONFIG" ]] && { echo "[ERROR] Config not found: $JOB_CONFIG"; exit 1; }

sbatch \
  --export=ALL \
  "$@" \
  "${GROMACS_SCRIPTS_DIR}/REMD-gromacs.sbatch" ${JOB_CONFIG:+"$JOB_CONFIG"}
