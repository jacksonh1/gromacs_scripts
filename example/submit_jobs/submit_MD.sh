#!/usr/bin/env bash
# Submit a plain production MD job. Copy this file, set your parameters, and run it.
# Usage: bash submit_MD.sh

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/gromacs_scripts"

# === Temperature & production ===
T_SIM=300        # single production temperature (K)
TOTAL_NS=2       # total production length (ns)
TRAJ_PS=10       # write a trajectory frame every this many ps

# === Equilibration (optional) ===
# EQ_NPT_NS: unrestrained NPT equilibration before production (ns).
#   0 (default) — restraints release at production START, so the trajectory
#                 captures the relaxation from the designed pose (stability /
#                 flexibility / variant comparison).
#   >0          — equilibrate first so production starts pre-relaxed (e.g. bound-
#                 state equilibrium sampling). Uncomment to enable:
# EQ_NPT_NS=0.5

# === System ===
PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"

OUTBASE="$(basename "${PDB_IN%.*}")"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_MD/"$OUTBASE"-"$TOTAL_NS"ns-MD-"$T_SIM"K-NPT"


sbatch \
  --export=ALL,PDB_IN="$PDB_IN",OUTBASE="$OUTBASE",OUTDIR="$OUTDIR",T_SIM="$T_SIM",TOTAL_NS="$TOTAL_NS",TRAJ_PS="$TRAJ_PS",EQ_NPT_NS="${EQ_NPT_NS:-0}" \
  "${GROMACS_SCRIPTS_DIR}/MD-gromacs.sbatch"
