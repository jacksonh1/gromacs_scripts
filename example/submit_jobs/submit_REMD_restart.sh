#!/usr/bin/env bash
# Resume a T-REMD production run that stopped early (TIMEOUT, node failure, scancel).
# Copy this file, set OUTDIR to the run you want to continue, and run it.
# Usage: bash submit_REMD_restart.sh

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

# The existing run directory to resume (contains prod/, logs/). Its
# prod/repNNN/remd.cpt checkpoints are what the restart continues from.
OUTDIR="/orcd/pool/004/jhalpin/09-fragfold/XBL_PEP_DESIGN/campaigns/2026-july-XBL_PEP_V2/target-TNFa/simulations/7kp7-TNFa-TNR1A_remd_959fba54"

# Match the ORIGINAL job's replica count (-n) and GPUs (--gres). Replica count
# is auto-detected from prod/rep* and asserted against -n inside the script.
REPLICAS=48

# Wall time for the *remaining* steps only. This run reached 16/20 ns in 8 h
# (~2 ns/h), so ~4 ns left ≈ 2 h — 4 h is a safe margin.
sbatch -n "$REPLICAS" --gres=gpu:l40s:4 -t 4:00:00 \
  --export=ALL,GROMACS_SCRIPTS_DIR="$GROMACS_SCRIPTS_DIR",OUTDIR="$OUTDIR" \
  "${GROMACS_SCRIPTS_DIR}/REMD-restart.sbatch"
