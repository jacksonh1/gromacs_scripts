#!/usr/bin/env bash
# Submit a REST2 (solute-tempering Hamiltonian replica exchange) job.
# Copy this file, set your parameters, and run it:  bash submit_REST2.sh
#
# REST2 runs every replica at the SAME physical temperature (T_MIN) but scales the
# PROTEIN force field by lambda_i = T_MIN/T_i across a geometric EFFECTIVE-temperature
# ladder (T_MIN..T_MAX). rep000 (lambda=1) is the true, unscaled T_MIN ensemble.
#
# Requires the GROMACS 2023.5 + PLUMED build (site_config.sh REST2_GMXRC) — the engine
# refuses to run without working hrex and fails loud if exchanges don't happen.

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

# === Replicas & effective-temperature ladder ===
REPLICAS=16
T_MIN=300          # physical temperature (ALL replicas)
T_MAX=400          # EFFECTIVE max solute temperature (top of the lambda ladder)

# === Production ===
TOTAL_NS=2
REPLEX_PS=1        # exchange attempt interval (ps); keep >= 1
ENSEMBLE=NPT       # NVT | NPT

# === System ===
PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"

OUTBASE="$(basename "${PDB_IN%.*}")"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_REST2/"$OUTBASE"-"$TOTAL_NS"ns-REST2-"$T_MIN"-"$T_MAX"Keff-"$REPLICAS"reps-"$ENSEMBLE"-exf-"$REPLEX_PS"ps"

# The engine's #SBATCH header targets pi_keating; add -p / --gres here to override,
# e.g. to build-portable mit_normal_gpu nodes:  -p mit_normal_gpu --gres=gpu:l40s:4
sbatch -n "$REPLICAS" --gres=gpu:2\
  --export=ALL,PDB_IN="$PDB_IN",OUTBASE="$OUTBASE",OUTDIR="$OUTDIR",T_MIN="$T_MIN",T_MAX="$T_MAX",TOTAL_NS="$TOTAL_NS",REPLEX_PS="$REPLEX_PS",ENSEMBLE="$ENSEMBLE" \
  "${GROMACS_SCRIPTS_DIR}/REST2-gromacs.sbatch"
