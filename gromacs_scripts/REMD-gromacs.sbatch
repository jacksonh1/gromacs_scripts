#!/bin/bash
#SBATCH --job-name=remd
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --gres=gpu:l40s:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# GROMACS T-REMD single-node workflow
# ---------------------------------------------------------------------------
# Intended usage:
#   Copy example/submit_jobs/submit_REMD.sh, set your parameters, run it.
#   That script sets GROMACS_SCRIPTS_DIR and passes parameters via --export.
#
# Advanced / direct usage:
#   sbatch --export=ALL,PDB_IN=my.pdb,... REMD-gromacs.sbatch
#   sbatch REMD-gromacs.sbatch config.sh   # source a config file as $1
###############################################################################

set -euo pipefail

CURRENT_TIME=$(date +%s)
CURRENT_STEP=0
PRESERVE_FROM_STEP=${PRESERVE_FROM_STEP:-9}

# ── Source optional config file ──────────────────────────────────────────────
if [[ -n "${1:-}" && -f "$1" ]]; then
  echo "[INFO] Sourcing config: $1"
  # shellcheck disable=SC1090
  source "$1"
fi

###############################################################################
# STEP 0 — Environment
###############################################################################
echo "===== JOB START $(date) ====="

: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — submit via example/submit_jobs/submit_REMD.sh}"
source "${GROMACS_SCRIPTS_DIR}/../site_config.sh"

module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"

set +u
source "$GMXRC"
set -u
GMX="${GMX:-gmx_mpi}"          # set GMX=gmx_mpi if your build uses that binary
MDRUN="${MDRUN:-gmx mdrun}"

# Verify GROMACS is available
command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }
echo "[OK] GROMACS found: $($GMX --version 2>&1 | head -1)"

###############################################################################
# STEP 1 — Parameters
###############################################################################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="

PDB_IN="${PDB_IN:-IL7-2-V3-cys.pdb}"
OUTBASE="${OUTBASE:-IL7-2-V3-cys}"
OUTDIR="${OUTDIR:-remd_${OUTBASE}}"
FORCE="${FORCE:-0}"

# Force field & water model
FF="${FF:-amber99sb-ildn}"          # GROMACS ff name
WATER="${WATER:-tip3p}"
BOX_SHAPE="${BOX_SHAPE:-dodecahedron}"   # dodecahedron ≈ truncated octahedron
BOX_BUFFER="${BOX_BUFFER:-1.0}"          # nm (GROMACS uses nm; 10 Å = 1.0 nm)
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

# Simulation parameters
REPLICAS=${SLURM_NTASKS:-${REPLICAS:-48}}
DT_PS="${DT_PS:-0.002}"
T_MIN="${T_MIN:-300}"
T_MAX="${T_MAX:-400}"
CUTOFF_NM="${CUTOFF_NM:-1.0}"           # 10 Å = 1.0 nm
GAMMA_LN="${GAMMA_LN:-2.0}"             # ps^-1  (tau_t = 1/gamma)
TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")

# REMD exchange & production
REPLEX_PS="${REPLEX_PS:-0.5}"            # attempt exchange every this many ps
TOTAL_NS="${TOTAL_NS:-20}"               # total production per replica (ns)
REPLEX_STEPS=$(python3 -c "print(int(round(float('${REPLEX_PS}')/float('${DT_PS}'))))")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")

# Equilibration
EQUI_NS="${EQUI_NS:-0.2}"               # NVT equil length per replica (ns)
EQUI_NSTEPS=$(python3 -c "print(int(round(float('${EQUI_NS}')*1000/float('${DT_PS}'))))")
NPT_SEG_STEPS="${NPT_SEG_STEPS:-10000}"
NPT_MIN_SEG="${NPT_MIN_SEG:-8}"
NPT_MAX_SEG="${NPT_MAX_SEG:-20}"
NPT_TOL_REL="${NPT_TOL_REL:-0.005}"

# Scratch — large trajectory files go here
SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
SCRATCH_DIR="${SCRATCH_DIR:-${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}}"
mkdir -p "$SCRATCH_DIR"

# Cleanup function
cleanup_scratch_on_error() {
  local rc=$?
  if [[ -n "$SCRATCH_DIR" && -d "$SCRATCH_DIR" ]]; then
    if (( CURRENT_STEP >= PRESERVE_FROM_STEP )); then
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); PRESERVING scratch: $SCRATCH_DIR"
    else
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); cleaning scratch: $SCRATCH_DIR"
      rm -rf "$SCRATCH_DIR"
    fi
  fi
  exit $rc
}
trap 'cleanup_scratch_on_error' ERR

# Check scratch space
SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if (( SCRATCH_AVAIL < 100 )); then
  echo "[ERROR] Low disk space in scratch: ${SCRATCH_AVAIL}GB available (need ~100GB)"
  exit 1
fi

# ── Output directory tree ────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"/{build,em,npt,nvt,prod,logs,trajectories} ./logs

[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set  (${REPLICAS} replicas, ${T_MIN}–${T_MAX} K, ${TOTAL_NS} ns)"

###############################################################################
# STEP 2 — Temperature ladder
###############################################################################
CURRENT_STEP=2
echo "===== STEP 2: Temperature ladder ====="

TEMPS=$(python3 - <<'PY'
import os, sys
n  = int(os.environ.get("REPLICAS", 48))
t0 = float(os.environ.get("T_MIN", 300))
t1 = float(os.environ.get("T_MAX", 400))
tl = os.environ.get("TEMPS_LIST", "").strip()
if tl:
    ts = [float(x) for x in tl.split(",")]
else:
    r = (t1/t0)**(1/(n-1)) if n > 1 else 1.0
    ts = [t0 * r**i for i in range(n)]
if len(ts) != n:
    print(f"[ERROR] {len(ts)} temps but {n} replicas", file=sys.stderr)
    sys.exit(1)
print(",".join(f"{t:.3f}" for t in ts))
PY
)
IFS=, read -r -a TARR <<< "$TEMPS"
echo "[OK] ${#TARR[@]} temperatures: ${TARR[0]} … ${TARR[-1]} K"

###############################################################################
# STEP 3 — Build system (pdb2gmx → editconf → solvate → genion)
###############################################################################
CURRENT_STEP=3
echo "===== STEP 3: Building system ====="

BUILD="${OUTDIR}/build"

# pdb2gmx — generate topology
echo "[INFO] Running pdb2gmx (ff=${FF}, water=${WATER})…"
$GMX pdb2gmx \
  -f "$PDB_IN" \
  -o "${BUILD}/${OUTBASE}_processed.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  -i "${BUILD}/${OUTBASE}_posre.itp" \
  -ff "$FF" -water "$WATER" -ignh \
  > "${BUILD}/pdb2gmx.log" 2>&1

# editconf — set box
echo "[INFO] Setting up ${BOX_SHAPE} box with ${BOX_BUFFER} nm buffer…"
$GMX editconf \
  -f "${BUILD}/${OUTBASE}_processed.gro" \
  -o "${BUILD}/${OUTBASE}_box.gro" \
  -bt "$BOX_SHAPE" -d "$BOX_BUFFER" -c \
  > "${BUILD}/editconf.log" 2>&1

# solvate
echo "[INFO] Solvating…"
$GMX solvate \
  -cp "${BUILD}/${OUTBASE}_box.gro" \
  -cs spc216.gro \
  -o "${BUILD}/${OUTBASE}_solv.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  > "${BUILD}/solvate.log" 2>&1

# genion — neutralize + add salt
# We need a dummy .tpr for genion, so write a minimal ions.mdp
cat > "${BUILD}/ions.mdp" <<EOF
integrator = steep
nsteps     = 0
EOF
$GMX grompp \
  -f "${BUILD}/ions.mdp" \
  -c "${BUILD}/${OUTBASE}_solv.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  -o "${BUILD}/ions.tpr" \
  -po "${BUILD}/mdout_ions.mdp" \
  -maxwarn 5 \
  > "${BUILD}/grompp_ions.log" 2>&1

GENION_FLAGS="-pname NA -nname CL"
if (( NEUTRALIZE )); then
  GENION_FLAGS+=" -neutral"
fi
if python3 -c "import sys; sys.exit(0 if float('${SALT_MOLAR}')>0 else 1)"; then
  GENION_FLAGS+=" -conc ${SALT_MOLAR}"
fi

echo "[INFO] Adding ions (neutralize=${NEUTRALIZE}, salt=${SALT_MOLAR} M)…"
echo SOL | $GMX genion \
  -s "${BUILD}/ions.tpr" \
  -o "${BUILD}/${OUTBASE}_ions.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  ${GENION_FLAGS} \
  > "${BUILD}/genion.log" 2>&1

# Final system files — these are our "topology" equivalents of prmtop/inpcrd
SYS_GRO="${BUILD}/${OUTBASE}_ions.gro"
SYS_TOP="${BUILD}/${OUTBASE}.top"
[[ -f "$SYS_GRO" && -f "$SYS_TOP" ]] || { echo "[ERROR] System build failed"; exit 1; }
echo "[OK] System built"

###############################################################################
# STEP 4 — Energy minimization
###############################################################################
CURRENT_STEP=4
echo "===== STEP 4: Energy minimization ====="

EM_DIR="${OUTDIR}/em"
cat > "${EM_DIR}/em.mdp" <<EOF
integrator  = steep
emtol       = 1000.0        ; kJ/mol/nm
emstep      = 0.01
nsteps      = 3000

nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
pbc         = xyz

; Position restraints on protein heavy atoms (uses posre.itp from pdb2gmx)
define      = -DPOSRES
EOF

$GMX grompp \
  -f "${EM_DIR}/em.mdp" \
  -c "$SYS_GRO" \
  -r "$SYS_GRO" \
  -p "$SYS_TOP" \
  -o "${EM_DIR}/em.tpr" \
  -po "${EM_DIR}/mdout_em.mdp" \
  -maxwarn 5 \
  > "${EM_DIR}/grompp_em.log" 2>&1

$GMX mdrun \
  -v -deffnm "${EM_DIR}/em" \
  -ntomp 8 \
  > "${EM_DIR}/mdrun_em.log" 2>&1

[[ -f "${EM_DIR}/em.gro" ]] || { echo "[ERROR] Minimization failed"; exit 1; }
echo "[OK] Minimization complete"

###############################################################################
# STEP 5 — NPT density equilibration (iterative convergence)
###############################################################################
CURRENT_STEP=5
echo "===== STEP 5: NPT density equilibration (iterative) ====="

NPT_DIR="${OUTDIR}/npt"
LAST_GRO="${EM_DIR}/em.gro"
prev_vol=""
seg=0

while (( seg < NPT_MAX_SEG )); do
  seg=$((seg+1))
  echo "[INFO] NPT segment $seg"

  if (( seg == 1 )); then
    GEN_VEL="yes"
    GEN_TEMP="${T_MIN}"
    CONT_FLAGS="gen-vel = yes
gen-temp = ${T_MIN}"
    SEG_DT="0.001"
    SEG_STEPS="${NPT_SEG_STEPS}"
  else
    CONT_FLAGS="gen-vel = no
continuation = yes"
    SEG_DT="${DT_PS}"
    SEG_STEPS="${NPT_SEG_STEPS}"
  fi

  cat > "${NPT_DIR}/npt_seg${seg}.mdp" <<EOF
integrator  = md
dt          = ${SEG_DT}
nsteps      = ${SEG_STEPS}
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_MIN} ${T_MIN}

pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 1.0
ref-p       = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

constraints = h-bonds
constraint-algorithm = LINCS

${CONT_FLAGS}

define      = -DPOSRES
EOF

  $GMX grompp \
    -f "${NPT_DIR}/npt_seg${seg}.mdp" \
    -c "$LAST_GRO" \
    -r "$LAST_GRO" \
    -p "$SYS_TOP" \
    -o "${NPT_DIR}/npt_seg${seg}.tpr" \
    -po "${NPT_DIR}/mdout_npt${seg}.mdp" \
    -maxwarn 5 \
    > "${NPT_DIR}/grompp_npt${seg}.log" 2>&1

  $GMX mdrun \
    -deffnm "${NPT_DIR}/npt_seg${seg}" \
    -x "${SCRATCH_DIR}/${OUTBASE}_npt_seg${seg}.xtc" \
    -ntomp 8 -nb gpu -pme gpu -bonded gpu \
    > "${NPT_DIR}/mdrun_npt${seg}.log" 2>&1

  # Extract volume
  echo Volume | $GMX energy \
    -f "${NPT_DIR}/npt_seg${seg}.edr" \
    -o "${NPT_DIR}/volume_seg${seg}.xvg" \
    > /dev/null 2>&1

  vol_avg=$(awk '!/^[#@]/ {sum+=$2; n++} END {if(n>0) printf "%.6f", sum/n}' \
            "${NPT_DIR}/volume_seg${seg}.xvg")
  echo "[INFO] Segment $seg average volume: ${vol_avg} nm^3"

  if [[ -n "$prev_vol" ]]; then
    rel_change=$(python3 -c "
prev=$prev_vol; curr=$vol_avg
print(abs(curr-prev)/prev if prev!=0 else 1.0)
")
    echo "[INFO] Relative volume change: $rel_change"

    if (( seg >= NPT_MIN_SEG )); then
      converged=$(python3 -c "print(0 if float('$rel_change')>float('${NPT_TOL_REL}') else 1)")
      if (( converged )); then
        echo "[INFO] Density converged at segment $seg"
        break
      fi
    fi
  fi

  prev_vol="$vol_avg"
  LAST_GRO="${NPT_DIR}/npt_seg${seg}.gro"
done

if (( seg == NPT_MAX_SEG )); then
  echo "[WARNING] Reached NPT_MAX_SEG=${NPT_MAX_SEG} without convergence; using last segment"
fi

NPT_GRO="${NPT_DIR}/npt_seg${seg}.gro"
NPT_CPT="${NPT_DIR}/npt_seg${seg}.cpt"
[[ -f "$NPT_GRO" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 6 — Prepare per-replica NVT equilibration directories
###############################################################################
CURRENT_STEP=6
echo "===== STEP 6: Preparing per-replica NVT equilibration ====="

NVT_DIR="${OUTDIR}/nvt"

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  T=${TARR[$i]}
  REPDIR="${NVT_DIR}/rep${idx}"
  mkdir -p "$REPDIR"

  cat > "${REPDIR}/nvt.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${EQUI_NSTEPS}
nstxout-compressed = 2500
nstlog      = 1000
nstenergy   = 1000

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T} ${T}

pcoupl      = no

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = yes
gen-temp    = ${T}

nstcomm     = 500
EOF

  $GMX grompp \
    -f "${REPDIR}/nvt.mdp" \
    -c "$NPT_GRO" \
    -r "$NPT_GRO" \
    -p "$SYS_TOP" \
    -o "${REPDIR}/nvt.tpr" \
    -po "${REPDIR}/mdout_nvt.mdp" \
    -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1

  # Pre-create symlink so deffnm writes the xtc directly to scratch
  ln -sf "${SCRATCH_DIR}/nvt_rep${idx}.xtc" "${REPDIR}/nvt.xtc"
done
echo "[OK] NVT .tpr files prepared for all replicas"

###############################################################################
# STEP 7 — Run per-replica NVT equilibration
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Running NVT equilibration ====="

# Build the -multidir list
MULTIDIR_NVT=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  MULTIDIR_NVT+="${NVT_DIR}/rep${idx} "
done

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir $MULTIDIR_NVT \
  -deffnm nvt \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_nvt_equil.log" 2>&1

# Verify all replicas finished
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${NVT_DIR}/rep${idx}/nvt.gro" ]] || { echo "[ERROR] NVT equil failed for replica $idx"; exit 1; }
done
echo "[OK] NVT equilibration done"

###############################################################################
# STEP 8 — Prepare REMD production inputs
###############################################################################
CURRENT_STEP=8
echo "===== STEP 8: Preparing REMD production inputs ====="

PROD_DIR="${OUTDIR}/prod"

# Output frequency — write trajectory every 10x exchange interval
NSTXOUT_COMP=$(python3 -c "print(int(${REPLEX_STEPS}*10))")

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  T=${TARR[$i]}
  REPDIR="${PROD_DIR}/rep${idx}"
  mkdir -p "$REPDIR"

  cat > "${REPDIR}/remd.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${PROD_NSTEPS}
nstxout-compressed = ${NSTXOUT_COMP}
nstlog      = ${NSTXOUT_COMP}
nstenergy   = ${NSTXOUT_COMP}

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T} ${T}

pcoupl      = no

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = no
continuation = yes

nstcomm     = 500
EOF

  $GMX grompp \
    -f "${REPDIR}/remd.mdp" \
    -c "${NVT_DIR}/rep${idx}/nvt.gro" \
    -t "${NVT_DIR}/rep${idx}/nvt.cpt" \
    -p "$SYS_TOP" \
    -o "${REPDIR}/remd.tpr" \
    -po "${REPDIR}/mdout_remd.mdp" \
    -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1

  # Pre-create symlink so deffnm writes the xtc directly to scratch
  ln -sf "${SCRATCH_DIR}/remd_rep${idx}.xtc" "${REPDIR}/remd.xtc"
done
echo "[OK] REMD .tpr files prepared"

###############################################################################
# STEP 9 — Run T-REMD production
###############################################################################
CURRENT_STEP=9
echo "===== STEP 9: Running T-REMD production ====="
echo "[INFO] ${REPLICAS} replicas × ${TOTAL_NS} ns, exchange every ${REPLEX_PS} ps"

MULTIDIR_PROD=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  MULTIDIR_PROD+="${PROD_DIR}/rep${idx} "
done

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir $MULTIDIR_PROD \
  -deffnm remd \
  -replex "$REPLEX_STEPS" \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_remd.log" 2>&1

REMD_EXIT=$?
if (( REMD_EXIT != 0 )); then
  echo "[ERROR] REMD production failed (exit=$REMD_EXIT)"
  exit 1
fi
echo "[OK] REMD production done"

###############################################################################
# STEP 10 — Finalize outputs & symlinks
###############################################################################
CURRENT_STEP=10
echo "===== STEP 10: Finalizing outputs ====="

# Symlink large trajectory files from scratch into the output tree
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  for prefix in nvt remd; do
    src="${SCRATCH_DIR}/${prefix}_rep${idx}.xtc"
    if [[ -f "$src" ]]; then
      ln -sf "$src" "${OUTDIR}/trajectories/${prefix}_rep${idx}.xtc"
    fi
  done
done

# Symlink NPT trajectories from scratch
for f in "${SCRATCH_DIR}/${OUTBASE}_npt_seg"*.xtc; do
  [[ -f "$f" ]] && ln -sf "$f" "${OUTDIR}/trajectories/$(basename "$f")"
done

# Export final PDB from replica 000
echo 0 | $GMX trjconv \
  -s "${PROD_DIR}/rep000/remd.tpr" \
  -f "${PROD_DIR}/rep000/remd.gro" \
  -o "${OUTDIR}/${OUTBASE}_final_rep000.pdb" \
  -pbc whole \
  > /dev/null 2>&1 || true

echo "[OK] Outputs finalized"

###############################################################################
# STEP 11 — Write parameters log
###############################################################################
CURRENT_STEP=11
ELAPSED_S=$(( $(date +%s) - CURRENT_TIME ))

for dest in "${OUTDIR}/parameters.txt" "${SCRATCH_DIR}/parameters.txt"; do
  cat > "$dest" <<EOF
Elapsed seconds: $ELAPSED_S
Job Name: ${SLURM_JOB_NAME:-remd}
Job ID: ${SLURM_JOB_ID:-N/A}
Replicas: $REPLICAS
Total ns per replica: $TOTAL_NS
Temperature range: ${T_MIN}–${T_MAX} K
Force field: $FF
Water model: $WATER
Box shape: $BOX_SHAPE
Box buffer (nm): $BOX_BUFFER
Neutralize: $NEUTRALIZE
Salt molarity: $SALT_MOLAR
Cutoff (nm): $CUTOFF_NM
tau_t (ps): $TAU_T
Time step (ps): $DT_PS
Exchange interval (ps): $REPLEX_PS
Exchange interval (steps): $REPLEX_STEPS
Production steps: $PROD_NSTEPS
Output directory: $OUTDIR
Scratch directory: $SCRATCH_DIR
GPUs: ${SLURM_GPUS_ON_NODE:-unknown}
NPT convergence segments: $seg
EOF
done

###############################################################################
# STEP 12 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:   $OUTDIR"
echo "[INFO] Trajectory links:   ${OUTDIR}/trajectories/"
echo "[INFO] Actual trajectories: $SCRATCH_DIR"
echo "[INFO] REMD exchange log:  ${PROD_DIR}/rep000/remd.log  (or check all rep dirs)"
echo ""
echo "Next steps:"
echo "  1. Check exchange acceptance rates:"
echo "       grep 'Repl  average probabilities' ${PROD_DIR}/rep000/remd.log"
echo ""
echo "  2. The 300 K trajectory is prod/rep000/remd.xtc — use directly for analysis."
echo "     (Each replica runs at a fixed temperature; coordinates exchange between them."
echo "      No demux step is needed to obtain the constant-temperature ensemble.)"
echo ""
echo "  3. Higher-T replicas: ${PROD_DIR}/rep001/remd.xtc, rep002/remd.xtc, etc."
echo ""
echo "  4. Analysis tools: see analysis_scripts/ (to be added)"
echo "=========================================="
