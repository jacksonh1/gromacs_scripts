#!/usr/bin/env bash
# Submit a T-REMD job. Copy this file, set your parameters, and run it.
# Usage: bash submit_REMD.sh

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

# === Replicas & temperatures ===
REPLICAS=48
T_MIN=300
T_MAX=450

# === Production ===
TOTAL_NS=2
REPLEX_PS=1
ENSEMBLE=NPT     # NVT (constant volume) | NPT (C-rescale pressure coupling)

# === Reproducibility ===
# -1 (default) = fresh random seed each run (each run an independent sample).
# Set to a non-negative integer to pin the RNG (velocity generation, V-rescale/C-rescale
# streams, replica exchange; per-replica seed = SEED+i). NOTE: on GPU this fixes the
# initial conditions and RNG streams, not a bit-for-bit identical trajectory.
SEED=-1

# === System ===
PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"

OUTBASE="$(basename "${PDB_IN%.*}")"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_T-REMD/"$OUTBASE"-"$TOTAL_NS"ns-REMD-"$T_MIN"-"$T_MAX"K-"$REPLICAS"reps-"$ENSEMBLE"-exf-"$REPLEX_PS"ps"


sbatch -n "$REPLICAS" \
  --export=ALL,PDB_IN="$PDB_IN",OUTBASE="$OUTBASE",OUTDIR="$OUTDIR",T_MIN="$T_MIN",T_MAX="$T_MAX",TOTAL_NS="$TOTAL_NS",REPLEX_PS="$REPLEX_PS",ENSEMBLE="$ENSEMBLE",SEED="$SEED" \
  "${GROMACS_SCRIPTS_DIR}/REMD-gromacs.sbatch"
