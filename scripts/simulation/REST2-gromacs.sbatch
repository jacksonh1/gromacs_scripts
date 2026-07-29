#!/bin/bash
#SBATCH --job-name=rest2
#SBATCH -p pi_keating
#SBATCH -C rocky8
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --gres=gpu:l40s:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# GROMACS REST2 (replica exchange with solute tempering) single-node workflow
# ---------------------------------------------------------------------------
# Solute tempering: every replica runs at the SAME physical temperature (T_MIN)
# but with the PROTEIN's force field scaled by lambda_i = T_MIN/T_i (a geometric
# EFFECTIVE-temperature ladder, T_MIN..T_MAX). rep000 has lambda=1 → the true,
# unscaled T_MIN ensemble (analyze it directly, like rep000 in T-REMD). Water and
# ions are unscaled. Exchange uses PLUMED Hamiltonian replica exchange (-hrex).
#
# REQUIRES a GROMACS build with WORKING hrex — the 2024.3-plumed build silently
# gives ZERO exchanges. This engine uses REST2_GMXRC (a GROMACS 2023.5 build) from
# site_config.sh, and runs a fail-loud pre-flight acceptance check before the real
# production so a wrong build can never masquerade as a working REST2 run.
#
# Intended usage:
#   Copy example/submit_jobs/submit_REST2.sh, set parameters, run it.
###############################################################################

set -euo pipefail

CURRENT_TIME=$(date +%s)
CURRENT_STEP=0
PRESERVE_SCRATCH_FROM="${PRESERVE_SCRATCH_FROM:-prod}"   # prod | density | always | never

# ── Source optional config file ──────────────────────────────────────────────
: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — submit via example/submit_jobs/submit_REST2.sh}"
VALIDATE_PARAMS="${GROMACS_SCRIPTS_DIR}/validate_params.py"
if [[ -n "${1:-}" ]]; then
  [[ -f "$1" ]] || { echo "[ERROR] config file not found: $1"; exit 1; }
  python3 "$VALIDATE_PARAMS" --engine rest2 --check-keys "$1"
  echo "[INFO] Sourcing config: $1"
  # shellcheck disable=SC1090
  source "$1"
fi
CFG_REPLICAS="${REPLICAS:-}"

###############################################################################
# STEP 0 — Environment (REST2: 2023.5 GROMACS + PLUMED + hcoll compat libs)
###############################################################################
echo "===== JOB START $(date) ====="
source "${GROMACS_SCRIPTS_DIR}/../../site_config.sh"

module purge
module load gcc
module load deprecated-modules          # provides REST2_CUDA_MODULE (cuda 12.4)
module load "$REST2_CUDA_MODULE"
module load "$OPENMPI_MODULE"

set +u
source "$REST2_GMXRC"                    # GROMACS 2023.5 build with working hrex
source "$PLUMED_SH"                      # PLUMED_KERNEL + LD_LIBRARY_PATH
set -u
# gmx_mpi AND the plumed CLI link libhcoll/libocoms; some GPU nodes lack /opt/mellanox.
export LD_LIBRARY_PATH="${HCOLL_COMPAT_DIR}:${LD_LIBRARY_PATH:-}"
GMX="${GMX:-gmx_mpi}"

command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }
# NB: use grep -c with pipefail OFF — `grep -q` closes the pipe on first match, mdrun -h
# then dies with SIGPIPE (rc 141), and `set -o pipefail` would turn that into a false negative.
set +o pipefail
_HREX_OK=$("$GMX" mdrun -h 2>&1 | grep -ci 'hrex')
set -o pipefail
(( _HREX_OK > 0 )) || { echo "[ERROR] this GROMACS has no -hrex (not PLUMED-patched for hrex)"; exit 1; }
command -v plumed >/dev/null 2>&1 || { echo "[ERROR] plumed CLI not found"; exit 1; }
plumed partial_tempering --help >/dev/null 2>&1 || plumed patch -l >/dev/null 2>&1 \
  || { echo "[ERROR] plumed CLI cannot run (missing libs? check HCOLL_COMPAT_DIR=$HCOLL_COMPAT_DIR)"; exit 1; }
echo "[OK] GROMACS: $($GMX --version 2>&1 | grep -m1 'GROMACS version')"
echo "[OK] PLUMED : $(plumed info --version 2>&1)"

###############################################################################
# STEP 1 — Parameters
###############################################################################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="

PDB_IN="${PDB_IN:-}"                          # required — no default
OUTBASE="${OUTBASE:-$(basename "${PDB_IN%.*}")}"
OUTDIR="${OUTDIR:-rest2_${OUTBASE}}"
FORCE="${FORCE:-0}"

FF="${FF:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
BOX_SHAPE="${BOX_SHAPE:-dodecahedron}"
BOX_BUFFER="${BOX_BUFFER:-1.0}"
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

REPLICAS=${SLURM_NTASKS:-${REPLICAS:-24}}

_NODE_CORES=$(( ${SLURM_NTASKS:-1} * ${SLURM_CPUS_PER_TASK:-1} ))
NTOMP_SERIAL="${NTOMP_SERIAL:-$(( _NODE_CORES < 8 ? _NODE_CORES : 8 ))}"

DT_PS="${DT_PS:-0.002}"
T_MIN="${T_MIN:-300}"                    # physical temperature (ALL replicas run here)
T_MAX="${T_MAX:-450}"                    # EFFECTIVE max solute temperature (top of the lambda ladder)
TEMPS_LIST="${TEMPS_LIST:-}"             # optional explicit effective-T ladder "300,308,..." (else geometric)
CUTOFF_NM="${CUTOFF_NM:-1.0}"
GAMMA_LN="${GAMMA_LN:-2.0}"

REPLEX_PS="${REPLEX_PS:-1.0}"            # exchange attempt interval (ps); keep >= 1 (GPU deadlock gotcha)
TOTAL_NS="${TOTAL_NS:-20}"

EQUIL_NS="${EQUIL_NS:-0.2}"
DENSITY_SEG_STEPS="${DENSITY_SEG_STEPS:-10000}"
DENSITY_MIN_SEG="${DENSITY_MIN_SEG:-8}"
DENSITY_MAX_SEG="${DENSITY_MAX_SEG:-20}"
DENSITY_TOL_REL="${DENSITY_TOL_REL:-0.005}"

ENSEMBLE="${ENSEMBLE:-NVT}"
REF_P="${REF_P:-1.0}"
TAU_P="${TAU_P:-1.0}"

export PDB_IN OUTBASE OUTDIR FORCE FF WATER BOX_SHAPE BOX_BUFFER NEUTRALIZE \
       SALT_MOLAR REPLICAS NTOMP_SERIAL DT_PS T_MIN T_MAX TEMPS_LIST CUTOFF_NM GAMMA_LN \
       REPLEX_PS TOTAL_NS EQUIL_NS DENSITY_SEG_STEPS DENSITY_MIN_SEG DENSITY_MAX_SEG \
       DENSITY_TOL_REL ENSEMBLE REF_P TAU_P PRESERVE_SCRATCH_FROM CFG_REPLICAS
python3 "$VALIDATE_PARAMS" --engine rest2

TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")
REPLEX_STEPS=$(python3 -c "print(int(round(float('${REPLEX_PS}')/float('${DT_PS}'))))")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")
EQUIL_NSTEPS=$(python3 -c "print(int(round(float('${EQUIL_NS}')*1000/float('${DT_PS}'))))")
PREFLIGHT_NSTEPS=$(( REPLEX_STEPS * 20 ))   # ~20 exchange attempts for the acceptance gate

case "$PRESERVE_SCRATCH_FROM" in
  always)  PRESERVE_FROM_STEP=0 ;;
  density) PRESERVE_FROM_STEP=5 ;;
  prod)    PRESERVE_FROM_STEP=11 ;;
  never)   PRESERVE_FROM_STEP=999 ;;
  *) echo "[ERROR] invalid PRESERVE_SCRATCH_FROM='$PRESERVE_SCRATCH_FROM'"; exit 1 ;;
esac

if [[ "$ENSEMBLE" == "NPT" ]]; then
  PCOUPL_BLOCK="pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = ${TAU_P}
ref-p       = ${REF_P}
compressibility = 4.5e-5"
else
  PCOUPL_BLOCK="pcoupl      = no"
fi

SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
SCRATCH_DIR="${SCRATCH_DIR:-${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}}"
mkdir -p "$SCRATCH_DIR"

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

SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if (( SCRATCH_AVAIL < 100 )); then
  echo "[ERROR] Low disk space in scratch: ${SCRATCH_AVAIL}GB available (need ~100GB)"; exit 1
fi

mkdir -p "${OUTDIR}"/{build,em,density,topol,equil,prod,logs,trajectories} ./logs

[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set  (${REPLICAS} replicas, ${T_MIN} K physical, solute ${T_MIN}–${T_MAX} K effective, ${TOTAL_NS} ns)"

###############################################################################
# STEP 2 — Effective-temperature ladder + lambda factors
###############################################################################
CURRENT_STEP=2
echo "===== STEP 2: Effective-temperature ladder + lambda ====="

LAMBDAS=$(python3 - <<'PY'
import os, sys
n  = int(os.environ.get("REPLICAS", 24))
t0 = float(os.environ.get("T_MIN", 300))
t1 = float(os.environ.get("T_MAX", 450))
tl = os.environ.get("TEMPS_LIST", "").strip()
if tl:
    ts = [float(x) for x in tl.split(",")]
else:
    r = (t1/t0)**(1/(n-1)) if n > 1 else 1.0
    ts = [t0 * r**i for i in range(n)]
if len(ts) != n:
    print(f"[ERROR] {len(ts)} temps but {n} replicas", file=sys.stderr); sys.exit(1)
# lambda_i = T_MIN / T_eff_i   (rep000: T_eff=T_MIN → lambda=1, unscaled)
print(",".join(f"{t0/t:.6f}" for t in ts))
PY
)
IFS=, read -r -a LARR <<< "$LAMBDAS"
echo "[OK] ${#LARR[@]} lambdas: ${LARR[0]} … ${LARR[-1]}  (T_eff ${T_MIN}…${T_MAX} K, all run at ${T_MIN} K)"

###############################################################################
# STEP 3 — Build system (pdb2gmx → editconf → solvate → genion)
###############################################################################
CURRENT_STEP=3
echo "===== STEP 3: Building system ====="
BUILD="${OUTDIR}/build"

echo "[INFO] Running pdb2gmx (ff=${FF}, water=${WATER})…"
$GMX pdb2gmx -f "$PDB_IN" -o "${BUILD}/${OUTBASE}_processed.gro" -p "${BUILD}/${OUTBASE}.top" \
  -i "${BUILD}/${OUTBASE}_posre.itp" -ff "$FF" -water "$WATER" -ignh > "${BUILD}/pdb2gmx.log" 2>&1

echo "[INFO] Setting up ${BOX_SHAPE} box with ${BOX_BUFFER} nm buffer…"
$GMX editconf -f "${BUILD}/${OUTBASE}_processed.gro" -o "${BUILD}/${OUTBASE}_box.gro" \
  -bt "$BOX_SHAPE" -d "$BOX_BUFFER" -c > "${BUILD}/editconf.log" 2>&1

echo "[INFO] Solvating…"
$GMX solvate -cp "${BUILD}/${OUTBASE}_box.gro" -cs spc216.gro \
  -o "${BUILD}/${OUTBASE}_solv.gro" -p "${BUILD}/${OUTBASE}.top" > "${BUILD}/solvate.log" 2>&1

cat > "${BUILD}/ions.mdp" <<EOF
integrator = steep
nsteps     = 0
EOF
$GMX grompp -f "${BUILD}/ions.mdp" -c "${BUILD}/${OUTBASE}_solv.gro" -p "${BUILD}/${OUTBASE}.top" \
  -o "${BUILD}/ions.tpr" -po "${BUILD}/mdout_ions.mdp" -maxwarn 5 > "${BUILD}/grompp_ions.log" 2>&1

GENION_FLAGS="-pname NA -nname CL"
(( NEUTRALIZE )) && GENION_FLAGS+=" -neutral"
if python3 -c "import sys; sys.exit(0 if float('${SALT_MOLAR}')>0 else 1)"; then
  GENION_FLAGS+=" -conc ${SALT_MOLAR}"
fi
echo "[INFO] Adding ions (neutralize=${NEUTRALIZE}, salt=${SALT_MOLAR} M)…"
echo SOL | $GMX genion -s "${BUILD}/ions.tpr" -o "${BUILD}/${OUTBASE}_ions.gro" \
  -p "${BUILD}/${OUTBASE}.top" ${GENION_FLAGS} > "${BUILD}/genion.log" 2>&1

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
emtol       = 1000.0
emstep      = 0.01
nsteps      = 3000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
pbc         = xyz
define      = -DPOSRES
EOF
$GMX grompp -f "${EM_DIR}/em.mdp" -c "$SYS_GRO" -r "$SYS_GRO" -p "$SYS_TOP" \
  -o "${EM_DIR}/em.tpr" -po "${EM_DIR}/mdout_em.mdp" -maxwarn 5 > "${EM_DIR}/grompp_em.log" 2>&1
$GMX mdrun -v -deffnm "${EM_DIR}/em" -ntomp "$NTOMP_SERIAL" > "${EM_DIR}/mdrun_em.log" 2>&1
[[ -f "${EM_DIR}/em.gro" ]] || { echo "[ERROR] Minimization failed"; exit 1; }
echo "[OK] Minimization complete"

###############################################################################
# STEP 5 — NPT density equilibration (iterative convergence)
###############################################################################
CURRENT_STEP=5
echo "===== STEP 5: NPT density equilibration (iterative) ====="
DENSITY_DIR="${OUTDIR}/density"
LAST_GRO="${EM_DIR}/em.gro"
prev_vol=""
seg=0
while (( seg < DENSITY_MAX_SEG )); do
  seg=$((seg+1))
  echo "[INFO] NPT segment $seg"
  if (( seg == 1 )); then
    CONT_FLAGS="gen-vel = yes
gen-temp = ${T_MIN}"
    SEG_DT="0.001"
  else
    CONT_FLAGS="gen-vel = no
continuation = yes"
    SEG_DT="${DT_PS}"
  fi
  cat > "${DENSITY_DIR}/density_seg${seg}.mdp" <<EOF
integrator  = md
dt          = ${SEG_DT}
nsteps      = ${DENSITY_SEG_STEPS}
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
  $GMX grompp -f "${DENSITY_DIR}/density_seg${seg}.mdp" -c "$LAST_GRO" -r "$LAST_GRO" -p "$SYS_TOP" \
    -o "${DENSITY_DIR}/density_seg${seg}.tpr" -po "${DENSITY_DIR}/mdout_density${seg}.mdp" -maxwarn 5 \
    > "${DENSITY_DIR}/grompp_density${seg}.log" 2>&1
  $GMX mdrun -deffnm "${DENSITY_DIR}/density_seg${seg}" \
    -x "${SCRATCH_DIR}/${OUTBASE}_density_seg${seg}.xtc" \
    -ntomp "$NTOMP_SERIAL" -nb gpu -pme gpu -bonded gpu > "${DENSITY_DIR}/mdrun_density${seg}.log" 2>&1
  echo Volume | $GMX energy -f "${DENSITY_DIR}/density_seg${seg}.edr" \
    -o "${DENSITY_DIR}/volume_seg${seg}.xvg" > /dev/null 2>&1
  vol_avg=$(awk '!/^[#@]/ {sum+=$2; n++} END {if(n>0) printf "%.6f", sum/n}' "${DENSITY_DIR}/volume_seg${seg}.xvg")
  echo "[INFO] Segment $seg average volume: ${vol_avg} nm^3"
  if [[ -n "$prev_vol" ]]; then
    rel_change=$(python3 -c "prev=$prev_vol; curr=$vol_avg; print(abs(curr-prev)/prev if prev!=0 else 1.0)")
    echo "[INFO] Relative volume change: $rel_change"
    if (( seg >= DENSITY_MIN_SEG )); then
      converged=$(python3 -c "print(0 if float('$rel_change')>float('${DENSITY_TOL_REL}') else 1)")
      if (( converged )); then echo "[INFO] Density converged at segment $seg"; break; fi
    fi
  fi
  prev_vol="$vol_avg"
  LAST_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
done
(( seg == DENSITY_MAX_SEG )) && echo "[WARNING] Reached DENSITY_MAX_SEG without convergence; using last segment"
DENSITY_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
[[ -f "$DENSITY_GRO" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 6 — Scaled topologies (grompp -pp → mark protein → partial_tempering)
###############################################################################
CURRENT_STEP=6
echo "===== STEP 6: Generating per-replica scaled topologies ====="
TOPOL_DIR="${OUTDIR}/topol"

# Fully-expanded processed topology (NO -DPOSRES → the scaled tops carry no restraints,
# matching equil/prod which are unrestrained).
cat > "${TOPOL_DIR}/pp.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = 0
cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
pbc         = xyz
EOF
$GMX grompp -f "${TOPOL_DIR}/pp.mdp" -c "$DENSITY_GRO" -p "$SYS_TOP" \
  -pp "${TOPOL_DIR}/processed.top" -o "${TOPOL_DIR}/pp.tpr" -po "${TOPOL_DIR}/mdout_pp.mdp" -maxwarn 5 \
  > "${TOPOL_DIR}/grompp_pp.log" 2>&1

echo "[INFO] Marking hot region (whole protein)…"
python3 "${GROMACS_SCRIPTS_DIR}/mark_hot_region.py" "${TOPOL_DIR}/processed.top" > "${TOPOL_DIR}/marked.top"

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  lam=${LARR[$i]}
  plumed partial_tempering "$lam" < "${TOPOL_DIR}/marked.top" > "${TOPOL_DIR}/topol_rep${idx}.top" 2> "${TOPOL_DIR}/pt_rep${idx}.log"
  [[ -s "${TOPOL_DIR}/topol_rep${idx}.top" ]] || { echo "[ERROR] partial_tempering failed for rep${idx} (lambda=$lam)"; exit 1; }
done
echo "[OK] ${REPLICAS} scaled topologies written (lambda ${LARR[0]}…${LARR[-1]})"

###############################################################################
# STEP 7 — Prepare per-replica equilibration (all at T_MIN, scaled Hamiltonians)
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Preparing per-replica equilibration (${ENSEMBLE}) ====="
EQUIL_DIR="${OUTDIR}/equil"
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  REPDIR="${EQUIL_DIR}/rep${idx}"; mkdir -p "$REPDIR"
  cat > "${REPDIR}/equil.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${EQUIL_NSTEPS}
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
ref-t       = ${T_MIN} ${T_MIN}
${PCOUPL_BLOCK}
constraints = h-bonds
constraint-algorithm = LINCS
gen-vel     = yes
gen-temp    = ${T_MIN}
nstcomm     = 500
EOF
  $GMX grompp -f "${REPDIR}/equil.mdp" -c "$DENSITY_GRO" -r "$DENSITY_GRO" \
    -p "${TOPOL_DIR}/topol_rep${idx}.top" -o "${REPDIR}/equil.tpr" -po "${REPDIR}/mdout_equil.mdp" -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  ln -sf "${SCRATCH_DIR}/equil_rep${idx}.xtc" "${REPDIR}/equil.xtc"
done
echo "[OK] equilibration .tpr files prepared"

###############################################################################
# STEP 8 — Run per-replica equilibration (independent MD, NO exchange)
###############################################################################
CURRENT_STEP=8
echo "===== STEP 8: Running per-replica equilibration ====="
MULTIDIR_EQUIL=""
for i in $(seq 0 $((REPLICAS-1))); do MULTIDIR_EQUIL+="${EQUIL_DIR}/rep$(printf %03d $i) "; done
mpirun -np "$REPLICAS" "$GMX" mdrun -multidir $MULTIDIR_EQUIL -deffnm equil \
  -ntomp 1 -nb gpu -pme gpu -bonded gpu -pin on > "${OUTDIR}/logs/mdrun_equil.log" 2>&1
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${EQUIL_DIR}/rep${idx}/equil.gro" ]] || { echo "[ERROR] equilibration failed for replica $idx"; exit 1; }
done
echo "[OK] per-replica equilibration done"

###############################################################################
# STEP 9 — Prepare REST2 production inputs
###############################################################################
CURRENT_STEP=9
echo "===== STEP 9: Preparing REST2 production inputs ====="
PROD_DIR="${OUTDIR}/prod"
NSTXOUT_COMP=$(python3 -c "print(int(${REPLEX_STEPS}*10))")
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  REPDIR="${PROD_DIR}/rep${idx}"; mkdir -p "$REPDIR"
  cat > "${REPDIR}/rest2.mdp" <<EOF
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
ref-t       = ${T_MIN} ${T_MIN}
${PCOUPL_BLOCK}
constraints = h-bonds
constraint-algorithm = LINCS
gen-vel     = no
continuation = yes
nstcomm     = 500
EOF
  $GMX grompp -f "${REPDIR}/rest2.mdp" -c "${EQUIL_DIR}/rep${idx}/equil.gro" -t "${EQUIL_DIR}/rep${idx}/equil.cpt" \
    -p "${TOPOL_DIR}/topol_rep${idx}.top" -o "${REPDIR}/rest2.tpr" -po "${REPDIR}/mdout_rest2.mdp" -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  printf "# REST2 hrex, no CVs (rep%s)\n" "$idx" > "${REPDIR}/plumed.dat"
  ln -sf "${SCRATCH_DIR}/rest2_rep${idx}.xtc" "${REPDIR}/rest2.xtc"
done
echo "[OK] REST2 .tpr files prepared"

###############################################################################
# STEP 10 — hrex pre-flight acceptance gate (FAIL LOUD if zero exchanges)
###############################################################################
CURRENT_STEP=10
echo "===== STEP 10: hrex pre-flight acceptance check (~20 exchange attempts) ====="
PRE_DIR="${SCRATCH_DIR}/preflight"
MULTIDIR_PRE=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  mkdir -p "${PRE_DIR}/rep${idx}"
  cp "${PROD_DIR}/rep${idx}/rest2.tpr" "${PROD_DIR}/rep${idx}/plumed.dat" "${PRE_DIR}/rep${idx}/"
  MULTIDIR_PRE+="${PRE_DIR}/rep${idx} "
done
mpirun -np "$REPLICAS" "$GMX" mdrun -multidir $MULTIDIR_PRE -deffnm rest2 \
  -plumed plumed.dat -hrex -replex "$REPLEX_STEPS" -nsteps "$PREFLIGHT_NSTEPS" \
  -ntomp 1 -nb gpu -pme gpu -bonded gpu -pin on -dlb no > "${OUTDIR}/logs/mdrun_preflight.log" 2>&1
echo "[CMD] python3 ${GROMACS_SCRIPTS_DIR}/check_hrex_acceptance.py ${PRE_DIR}/rep000/rest2.log"
python3 "${GROMACS_SCRIPTS_DIR}/check_hrex_acceptance.py" "${PRE_DIR}/rep000/rest2.log"
rm -rf "$PRE_DIR"
echo "[OK] hrex is exchanging — proceeding to production"

###############################################################################
# STEP 11 — Run REST2 production
###############################################################################
CURRENT_STEP=11
echo "===== STEP 11: Running REST2 production ====="
echo "[INFO] ${REPLICAS} replicas × ${TOTAL_NS} ns, exchange every ${REPLEX_PS} ps, ${T_MIN} K physical"
MULTIDIR_PROD=""
for i in $(seq 0 $((REPLICAS-1))); do MULTIDIR_PROD+="${PROD_DIR}/rep$(printf %03d $i) "; done
mpirun -np "$REPLICAS" "$GMX" mdrun -multidir $MULTIDIR_PROD -deffnm rest2 \
  -plumed plumed.dat -hrex -replex "$REPLEX_STEPS" \
  -ntomp 1 -nb gpu -pme gpu -bonded gpu -pin on -dlb no > "${OUTDIR}/logs/mdrun_rest2.log" 2>&1
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${PROD_DIR}/rep${idx}/rest2.gro" ]] || { echo "[ERROR] REST2 production failed for replica $idx"; exit 1; }
done
echo "[OK] REST2 production done"

###############################################################################
# STEP 12 — Finalize outputs & symlinks
###############################################################################
CURRENT_STEP=12
echo "===== STEP 12: Finalizing outputs ====="
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  for prefix in equil rest2; do
    src="${SCRATCH_DIR}/${prefix}_rep${idx}.xtc"
    [[ -f "$src" ]] && ln -sf "$src" "${OUTDIR}/trajectories/${prefix}_rep${idx}.xtc"
  done
done
for f in "${SCRATCH_DIR}/${OUTBASE}_density_seg"*.xtc; do
  [[ -f "$f" ]] && ln -sf "$f" "${OUTDIR}/trajectories/$(basename "$f")"
done
echo 0 | $GMX trjconv -s "${PROD_DIR}/rep000/rest2.tpr" -f "${PROD_DIR}/rep000/rest2.gro" \
  -o "${OUTDIR}/${OUTBASE}_final_rep000.pdb" -pbc whole > /dev/null 2>&1 || true
echo "[OK] Outputs finalized"

###############################################################################
# STEP 13 — Write parameters log
###############################################################################
CURRENT_STEP=13
ELAPSED_S=$(( $(date +%s) - CURRENT_TIME ))
for dest in "${OUTDIR}/parameters.txt" "${SCRATCH_DIR}/parameters.txt"; do
  cat > "$dest" <<EOF
Elapsed seconds: $ELAPSED_S
Job Name: ${SLURM_JOB_NAME:-rest2}
Job ID: ${SLURM_JOB_ID:-N/A}
Method: REST2 (solute tempering, PLUMED hrex)
Replicas: $REPLICAS
Physical temperature (all replicas): ${T_MIN} K
Effective solute-temperature ladder: ${T_MIN}–${T_MAX} K
Lambda ladder: ${LARR[0]} … ${LARR[-1]}
Hot region: whole protein
Total ns per replica: $TOTAL_NS
Force field: $FF
Water model: $WATER
Box shape: $BOX_SHAPE ; buffer (nm): $BOX_BUFFER
Neutralize: $NEUTRALIZE ; Salt molarity: $SALT_MOLAR
Cutoff (nm): $CUTOFF_NM ; tau_t (ps): $TAU_T ; Time step (ps): $DT_PS
Exchange interval: ${REPLEX_PS} ps (${REPLEX_STEPS} steps)
Production ensemble: $ENSEMBLE$( [[ "$ENSEMBLE" == NPT ]] && echo " (C-rescale, ref-p ${REF_P} bar, tau-p ${TAU_P} ps)" )
Per-replica equil ns: $EQUIL_NS
GROMACS: $($GMX --version 2>&1 | grep -m1 'GROMACS version')
PLUMED: $(plumed info --version 2>&1)
Output directory: $OUTDIR
Scratch directory: $SCRATCH_DIR
GPUs: ${SLURM_GPUS_ON_NODE:-unknown}
Density convergence segments: $seg
EOF
done

###############################################################################
# STEP 14 — Post-analysis (rep000 = lambda=1 = physical T_MIN ensemble)
###############################################################################
CURRENT_STEP=14
echo "===== STEP 14: Post-analysis (rep000 / lambda=1 / ${T_MIN} K) ====="
ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR 000"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" 000 \
  || echo "[WARN] run_analysis.sh failed — re-run the command above"

###############################################################################
# STEP 15 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:   $OUTDIR"
echo "[INFO] Trajectory links:   ${OUTDIR}/trajectories/"
echo "[INFO] Actual trajectories: $SCRATCH_DIR"
echo "[INFO] rep000 = lambda=1 = ${T_MIN} K physical ensemble (analyze directly)"
echo "[INFO] Exchange log:       ${PROD_DIR}/rep000/rest2.log"
echo "[INFO] Acceptance rates:   ${OUTDIR}/analysis/remd_acceptance.csv"
echo "=========================================="
