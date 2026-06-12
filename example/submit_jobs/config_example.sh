#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Example config file for REMD-gromacs.sbatch
# ─────────────────────────────────────────────────────────────────────────────
# Usage:  sbatch REMD-gromacs.sbatch config.sh
#
# Any variable left commented here will use the default in the main script.
# ─────────────────────────────────────────────────────────────────────────────

# ── Replicas & temperatures ─────────────────────────────────────────────────
REPLICAS=48
T_MIN=300
T_MAX=450
TOTAL_NS=2                           # ns per replica
REPLEX_PS=1                       # exchange attempt interval (ps)
ENSEMBLE=NVT                         # NVT (constant volume) | NPT (C-rescale barostat)
# REF_P="1.0"                         # bar (NPT only)
# TAU_P="1.0"                         # ps  (NPT only)
# TEMPS_LIST=""                       # override: "300.0,305.2,310.5,..."

# ── Input ────────────────────────────────────────────────────────────────────
PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"
OUTBASE="$(basename "${PDB_IN%.*}")"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_T-REMD/"$OUTBASE"-"$TOTAL_NS"ns-REMD-"$T_MIN"-"$T_MAX"K-"$REPLICAS"reps-"$ENSEMBLE"-exf-"$REPLEX_PS"ps"
# OUTDIR="remd_${OUTBASE}"            # auto-derived from OUTBASE

# ── Force field / box ───────────────────────────────────────────────────────
# FF="amber14sb"                      # pdb2gmx force field name
# WATER="tip3p"
# BOX_SHAPE="dodecahedron"            # dodecahedron ≈ truncated octahedron
# BOX_BUFFER="1.0"                    # nm  (10 Å)
# NEUTRALIZE=1
# SALT_MOLAR="0.15"


# ── MD parameters ───────────────────────────────────────────────────────────
# DT_PS="0.002"                       # timestep in ps
# CUTOFF_NM="0.9"                     # non-bonded cutoff in nm (9 Å)
# GAMMA_LN="2.0"                      # Langevin friction (ps^-1)


# ── Equilibration ────────────────────────────────────────────────────────────
# EQUIL_NS="0.2"                      # per-replica equilibration (ns)
# DENSITY_SEG_STEPS=10000             # steps per density-equilibration segment
# DENSITY_MIN_SEG=8                   # min segments before convergence check
# DENSITY_MAX_SEG=20
# DENSITY_TOL_REL="0.005"             # relative volume change tolerance

# ── Scratch ──────────────────────────────────────────────────────────────────
# SCRATCH_DIR="/path/to/fast/storage/${SLURM_JOB_ID}"
# PRESERVE_SCRATCH_FROM=prod          # prod|density|always|never — keep scratch on failure from this stage

# ── GROMACS binary ───────────────────────────────────────────────────────────
# GMX="gmx_mpi"                       # this build provides gmx_mpi only; leave unset
