#!/usr/bin/env bash
# Submit a plain production MD job. Copy this file, set your parameters, and run it.
# Usage: bash submit_MD.sh

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

# === Temperature & production ===
T_SIM=300        # single production temperature (K)
TOTAL_NS=2       # total production length (ns)
TRAJ_PS=10       # write a trajectory frame every this many ps

# === Reproducibility ===
# -1 (default) = fresh random seed each run. Set to a non-negative integer to pin the RNG
# (velocity generation + V-rescale/C-rescale streams). NOTE: on GPU this fixes the initial
# conditions and RNG streams, not a bit-for-bit identical trajectory.
SEED=-1

# === Equilibration (optional) ===
# RELAX_NS: unrestrained NPT relaxation (relax/ stage) before production (ns).
#   0 (default) — restraints release at production START, so the trajectory
#                 captures the relaxation from the designed pose (stability /
#                 flexibility / variant comparison).
#   >0          — equilibrate first so production starts pre-relaxed (e.g. bound-
#                 state equilibrium sampling). Uncomment to enable:
# RELAX_NS=0.5

# === Force field ===
# Default is AMBER99SB-ILDN. For CHARMM36m use the alias "charmm36m" — the engine
# resolves it to the installed dated port (charmm36-feb2026_cgenff-5.0) and records
# that resolved name in parameters.txt. Aliases are defined in site_config.sh.
# CHARMM auto-switches the mdp to force-switch vdW at 1.2 nm and forces CUTOFF_NM=1.2.
FF="amber99sb-ildn"     # or: FF="charmm36m"
WATER="tip3p"           # resolved inside the FF dir; always matches the force field

# === System ===
PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"

OUTBASE="$(basename "${PDB_IN%.*}")"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_MD/"$OUTBASE"-"$TOTAL_NS"ns-MD-"$T_SIM"K-NPT"


sbatch \
  --export=ALL,PDB_IN="$PDB_IN",OUTBASE="$OUTBASE",OUTDIR="$OUTDIR",T_SIM="$T_SIM",TOTAL_NS="$TOTAL_NS",TRAJ_PS="$TRAJ_PS",RELAX_NS="${RELAX_NS:-0}",SEED="$SEED",FF="$FF",WATER="$WATER" \
  "${GROMACS_SCRIPTS_DIR}/MD-gromacs.sbatch"
