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
# FF="amber14sb"                      # pdb2gmx force field name (default: amber99sb-ildn)
#                                     #   FF="charmm36m" → CHARMM36m, via the alias defined
#                                     #   in site_config.sh (FF_ALIASES). The engine resolves
#                                     #   it to the installed dated port name and logs THAT in
#                                     #   parameters.txt, so the release is always on record.
#                                     #   The real directory name works too. Any FF starting
#                                     #   "charmm" auto-switches the mdp to force-switch vdW
#                                     #   at 1.2 nm + DispCorr=no and forces CUTOFF_NM=1.2.
# WATER="tip3p"                       # resolved inside the FF dir: amber→standard TIP3P,
#                                     #   charmm→CHARMM-modified TIP3P (automatic, matched)
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
# SEED=-1                             # -1 (default): random per run. >=0: pin RNG (velocity gen,
#                                     #   V-rescale/C-rescale, exchange; per-replica = SEED+i).
#                                     #   GPU runs are not bit-for-bit reproducible regardless.

# ── Production ───────────────────────────────────────────────────────────────
TOTAL_NS=20                           # ns per replica
REPLEX_PS="1.0"                       # exchange attempt interval (ps). Keep >= 1.0:
                                      #   sub-ps exchange deadlocks GPU-resident REMD
                                      #   (MPI collective hang, clean physics). See GOTCHAS.md.
# ENSEMBLE=NVT                        # NVT (constant volume) | NPT (C-rescale barostat)
# REF_P="1.0"                         # bar (NPT only)
# TAU_P="1.0"                         # ps  (NPT only)

# ── Equilibration ────────────────────────────────────────────────────────────
# EQUIL_NS="0.2"                      # per-replica equilibration (ns)
# DENSITY_SEG_STEPS=10000             # steps per density-equilibration segment
# DENSITY_MIN_SEG=8                   # min segments before convergence check
# DENSITY_MAX_SEG=20
# DENSITY_TOL_REL="0.005"             # relative volume change tolerance

# ── Scratch / output model ────────────────────────────────────────────────────
# SCRATCH_DIR="/path/to/fast/storage/${SLURM_JOB_ID}"
# PRESERVE_SCRATCH_FROM=prod          # prod|density|always|never — keep scratch on failure from this stage
# SYMLINK_BULK=1                      # 1 (default): bulk stage dirs (prod/ equil/ density/ …) are folder
#                                     #   symlinks into scratch — quota-safe, laptop-sync-friendly.
#                                     # 0: everything real in OUTDIR, no scratch offload (use when OUTDIR
#                                     #   is already on a large disk).

# ── GROMACS binary ───────────────────────────────────────────────────────────
# GMX="gmx_mpi"                       # this build provides gmx_mpi only; leave unset
