#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Example config file for REMD-gromacs.sbatch
# ─────────────────────────────────────────────────────────────────────────────
# Usage:  sbatch REMD-gromacs.sbatch config.sh
#
# Any variable left commented here will use the default in the main script.
# ─────────────────────────────────────────────────────────────────────────────

# ── Input ────────────────────────────────────────────────────────────────────
PDB_IN="IL7-2-V3-cys.pdb"
OUTBASE="IL7-2-V3-cys"
# OUTDIR="remd_${OUTBASE}"            # auto-derived from OUTBASE

# ── Force field / box ───────────────────────────────────────────────────────
# FF="amber14sb"                      # pdb2gmx force field name
# WATER="tip3p"
# BOX_SHAPE="dodecahedron"            # dodecahedron ≈ truncated octahedron
# BOX_BUFFER="1.0"                    # nm  (10 Å)
# NEUTRALIZE=1
# SALT_MOLAR="0.15"

# ── Replicas & temperatures ─────────────────────────────────────────────────
REPLICAS=48
T_MIN=300
T_MAX=400
# TEMPS_LIST=""                       # override: "300.0,305.2,310.5,..."

# ── MD parameters ───────────────────────────────────────────────────────────
# DT_PS="0.002"                       # timestep in ps
# CUTOFF_NM="0.9"                     # non-bonded cutoff in nm (9 Å)
# GAMMA_LN="2.0"                      # Langevin friction (ps^-1)

# ── Production ───────────────────────────────────────────────────────────────
TOTAL_NS=20                           # ns per replica
REPLEX_PS="0.5"                       # exchange attempt interval (ps)

# ── Equilibration ────────────────────────────────────────────────────────────
# EQUI_NS="0.2"                       # NVT equil per replica (ns)
# NPT_SEG_STEPS=10000                 # steps per NPT density segment
# NPT_MIN_SEG=8                       # min segments before convergence check
# NPT_MAX_SEG=20
# NPT_TOL_REL="0.005"                 # relative volume change tolerance

# ── Scratch ──────────────────────────────────────────────────────────────────
# SCRATCH_DIR="/path/to/fast/storage/${SLURM_JOB_ID}"
# PRESERVE_FROM_STEP=9                # preserve scratch on failure from this step

# ── GROMACS binary ───────────────────────────────────────────────────────────
# GMX="gmx"                           # or "gmx_mpi" depending on your build
# MDRUN="gmx mdrun"                   # the mdrun command (script appends _mpi)
