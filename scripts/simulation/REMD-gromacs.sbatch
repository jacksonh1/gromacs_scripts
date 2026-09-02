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
PRESERVE_SCRATCH_FROM="${PRESERVE_SCRATCH_FROM:-prod}"   # prod | density | always | never (resolved to a step index in STEP 1)

# ── Source optional config file ──────────────────────────────────────────────
: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — submit via example/submit_jobs/submit_REMD.sh}"
VALIDATE_PARAMS="${GROMACS_SCRIPTS_DIR}/validate_params.py"
if [[ -n "${1:-}" ]]; then
  [[ -f "$1" ]] || { echo "[ERROR] config file not found: $1"; exit 1; }
  python3 "$VALIDATE_PARAMS" --engine remd --check-keys "$1"   # reject typo'd keys before sourcing
  echo "[INFO] Sourcing config: $1"
  # shellcheck disable=SC1090
  source "$1"
fi
CFG_REPLICAS="${REPLICAS:-}"   # config/--export REPLICAS, before SLURM -n (ntasks) overrides it below

###############################################################################
# STEP 0 — Environment
###############################################################################
echo "===== JOB START $(date) ====="

source "${GROMACS_SCRIPTS_DIR}/../../site_config.sh"

module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"

set +u
source "$GMXRC"
set -u
GMX="${GMX:-gmx_mpi}"          # this build provides gmx_mpi only

# Verify GROMACS is available
command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }
echo "[OK] GROMACS found: $($GMX --version 2>&1 | head -1)"

###############################################################################
# STEP 1 — Parameters
###############################################################################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="

PDB_IN="${PDB_IN:-}"                          # required (validated below) — no default
OUTBASE="${OUTBASE:-$(basename "${PDB_IN%.*}")}"
OUTDIR="${OUTDIR:-remd_${OUTBASE}}"
FORCE="${FORCE:-0}"
# Output model: 1 (default) = bulk stage dirs (density/ equil/ prod/) are folder
# symlinks into scratch (quota-safe, laptop-sync-friendly); 0 = everything real in
# OUTDIR with no scratch offload.
SYMLINK_BULK="${SYMLINK_BULK:-1}"

# Force field & water model
FF="${FF:-amber99sb-ildn}"          # GROMACS ff name
# Resolve a force-field alias (site_config.sh FF_ALIASES) to the real directory
# name, e.g. charmm36m -> charmm36-feb2026_cgenff-5.0. A non-alias passes through.
# FF_REQUESTED keeps what the job asked for so parameters.txt can show both.
FF_REQUESTED="$FF"
FF="$(resolve_ff "$FF")"
if [[ $FF != "$FF_REQUESTED" ]]; then
  FF_ALIAS_NOTE="   (alias: ${FF_REQUESTED})"
  echo "[INFO] FF alias '${FF_REQUESTED}' resolved to force field '${FF}'"
else
  FF_ALIAS_NOTE=""
fi
WATER="${WATER:-tip3p}"
BOX_SHAPE="${BOX_SHAPE:-dodecahedron}"   # dodecahedron ≈ truncated octahedron
BOX_BUFFER="${BOX_BUFFER:-1.0}"          # nm (GROMACS uses nm; 10 Å = 1.0 nm)
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

# Simulation parameters
REPLICAS=${SLURM_NTASKS:-${REPLICAS:-48}}

# OpenMP threads for the SERIAL prep steps (EM, NPT density). These run as a
# single rank before the replicas exist, so they may use the whole node: derive
# the count from the full allocation (ntasks × cpus-per-task), capped at 8 to stay
# in GROMACS' good OpenMP-scaling range and to never request more threads than
# cores. Overridable via NTOMP_SERIAL. (The multi-replica equilibration/production steps use
# one rank per replica with -ntomp 1 instead — see STEPS 7 and 9.)
_NODE_CORES=$(( ${SLURM_NTASKS:-1} * ${SLURM_CPUS_PER_TASK:-1} ))
NTOMP_SERIAL="${NTOMP_SERIAL:-$(( _NODE_CORES < 8 ? _NODE_CORES : 8 ))}"

DT_PS="${DT_PS:-0.002}"
T_MIN="${T_MIN:-300}"
T_MAX="${T_MAX:-400}"
CUTOFF_NM="${CUTOFF_NM:-1.0}"           # 10 Å = 1.0 nm

# ── Force-field-dependent nonbonded settings ─────────────────────────────────
# CHARMM force fields require force-switched vdW (switch at 1.0 nm, rvdw=rcoulomb=
# 1.2 nm) and NO analytic dispersion correction; AMBER uses a plain Verlet cutoff
# at CUTOFF_NM with DispCorr=EnerPres. FF=charmm* selects the CHARMM nonbonded
# block. VDW_BLOCK carries its own trailing newline, so the AMBER branch (empty)
# leaves the generated mdp byte-identical to before.
if [[ "$FF" == charmm* ]]; then
  [[ "$CUTOFF_NM" == "1.2" ]] || echo "[INFO] FF=$FF requires force-switched vdW at 1.2 nm — overriding CUTOFF_NM=${CUTOFF_NM} → 1.2 (switch at 1.0)"
  CUTOFF_NM="1.2"
  DISPCORR="no"
  VDW_BLOCK=$'vdwtype      = cutoff\nvdw-modifier = force-switch\nrvdw-switch  = 1.0\n'
else
  DISPCORR="EnerPres"
  VDW_BLOCK=""
fi
GAMMA_LN="${GAMMA_LN:-2.0}"             # ps^-1  (tau_t = 1/gamma)

# RNG seed. -1 (default) = GROMACS draws a fresh random seed per run (legacy behaviour,
# each run an independent sample). Set SEED>=0 to pin the setup: velocity generation,
# the stochastic V-rescale/C-rescale streams, and the replica-exchange RNG. Per-replica
# stages offset it by the replica index (SEED+i) so replicas stay decorrelated. NOTE: GPU
# runs are still not bit-for-bit reproducible (nondeterministic float reduction order); a
# fixed seed pins the initial conditions and RNG streams, not the exact trajectory.
SEED="${SEED:--1}"
# Per-replica seed: keep -1 (random) when SEED is -1, else offset the base by the replica
# index so each replica's velocity draw and V-rescale/C-rescale stream is decorrelated.
rep_seed() { if [[ "$SEED" == "-1" ]]; then echo -1; else echo $(( SEED + $1 )); fi; }

# REMD exchange & production
REPLEX_PS="${REPLEX_PS:-1.0}"            # exchange attempt interval (ps); keep >= 1 (GPU deadlock gotcha)
TOTAL_NS="${TOTAL_NS:-20}"               # total production per replica (ns)

# Equilibration
EQUIL_NS="${EQUIL_NS:-0.2}"               # per-replica equilibration length (ns)
DENSITY_SEG_STEPS="${DENSITY_SEG_STEPS:-10000}"
DENSITY_MIN_SEG="${DENSITY_MIN_SEG:-8}"
DENSITY_MAX_SEG="${DENSITY_MAX_SEG:-20}"
DENSITY_TOL_REL="${DENSITY_TOL_REL:-0.005}"

# Ensemble for production + per-replica equilibration
ENSEMBLE="${ENSEMBLE:-NVT}"              # NVT (constant volume) | NPT (C-rescale barostat)
REF_P="${REF_P:-1.0}"                    # bar  (NPT only)
TAU_P="${TAU_P:-1.0}"                    # ps   (NPT only)

# Validate resolved parameters before deriving anything or building. Export the job
# params so validate_params.py reads them from the environment (covers config- and
# --export-set values alike). Keep this list in sync with validate_params.py.
export PDB_IN OUTBASE OUTDIR FORCE FF WATER BOX_SHAPE BOX_BUFFER NEUTRALIZE \
       SALT_MOLAR REPLICAS NTOMP_SERIAL DT_PS T_MIN T_MAX CUTOFF_NM GAMMA_LN \
       REPLEX_PS TOTAL_NS EQUIL_NS DENSITY_SEG_STEPS DENSITY_MIN_SEG DENSITY_MAX_SEG \
       DENSITY_TOL_REL ENSEMBLE REF_P TAU_P PRESERVE_SCRATCH_FROM SYMLINK_BULK SEED CFG_REPLICAS
python3 "$VALIDATE_PARAMS" --engine remd

# Derived step counts (after validation, so a non-numeric input fails with a clear
# message above instead of a raw python traceback here).
TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")
REPLEX_STEPS=$(python3 -c "print(int(round(float('${REPLEX_PS}')/float('${DT_PS}'))))")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")
EQUIL_NSTEPS=$(python3 -c "print(int(round(float('${EQUIL_NS}')*1000/float('${DT_PS}'))))")

# Resolve the scratch-preservation stage name to this engine's step index (stage names
# mean the same thing across MD/REMD; the index they map to is engine-specific). The
# validator already restricted PRESERVE_SCRATCH_FROM to {prod,density,always,never}.
case "$PRESERVE_SCRATCH_FROM" in
  always)  PRESERVE_FROM_STEP=0 ;;
  density) PRESERVE_FROM_STEP=5 ;;
  prod)    PRESERVE_FROM_STEP=9 ;;
  never)   PRESERVE_FROM_STEP=999 ;;
  *) echo "[ERROR] invalid PRESERVE_SCRATCH_FROM='$PRESERVE_SCRATCH_FROM'"; exit 1 ;;
esac

# Pressure-coupling block injected into the per-replica equilibration (STEP 6) and
# production (STEP 8) mdp. NPT uses C-rescale — robust + correct ensemble, matching the
# density-equilibration barostat. GROMACS adds the PV term to the replica-exchange
# Metropolis criterion automatically when pressure coupling is on, so -replex is unchanged.
if [[ "$ENSEMBLE" == "NPT" ]]; then
  PCOUPL_BLOCK="pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = ${TAU_P}
ref-p       = ${REF_P}
compressibility = 4.5e-5"
else
  PCOUPL_BLOCK="pcoupl      = no"
fi

# Scratch — large trajectory files go here (only when SYMLINK_BULK=1). When 0, bulk
# stays real in OUTDIR: SCRATCH_DIR is left empty and every scratch-aware step below
# (ERR-trap cleanup, df check, parameters dual-write, finalize) skips it.
if (( SYMLINK_BULK )); then
  SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
  SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
  SCRATCH_DIR="${SCRATCH_DIR:-${SCRATCH_ROOT}/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}}"
  mkdir -p "$SCRATCH_DIR"
  DF_TARGET="$SCRATCH_DIR"
else
  SCRATCH_DIR=""                 # no scratch offload; bulk lives in OUTDIR
  mkdir -p "$OUTDIR"             # so the df check below has a real path to probe
  DF_TARGET="$OUTDIR"
fi

# Cleanup function — a no-op when SCRATCH_DIR is empty (SYMLINK_BULK=0: bulk is real
# in OUTDIR, which is never auto-deleted, so a failed run is preserved regardless).
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

# Check disk space where the bulk will actually land (scratch, or OUTDIR when SYMLINK_BULK=0)
SCRATCH_AVAIL=$(df -BG "$DF_TARGET" | awk 'NR==2 {print $4}' | sed 's/G//')
if (( SCRATCH_AVAIL < 100 )); then
  echo "[ERROR] Low disk space in ${DF_TARGET}: ${SCRATCH_AVAIL}GB available (need ~100GB)"
  exit 1
fi

# ── Output directory tree (folder-symlink model) ─────────────────────────────
# Small, laptop-worthy dirs stay REAL in OUTDIR. The bulk stage dirs (density/
# equil/ prod/) are SYMLINKS into scratch when SYMLINK_BULK=1 — created BEFORE the
# run so mdrun writes trajectories/intermediates straight onto scratch and the
# tight-quota pool is never touched — or REAL dirs in OUTDIR when SYMLINK_BULK=0.
# Either way the stage-dir variables below point at ${OUTDIR}/<stage>, so nothing
# downstream (deffnm/-x/grompp) changes. build/ stays real so the solvate/genion
# topology rename() stays same-filesystem (see gotcha).
mkdir -p "${OUTDIR}"/{build,em,logs,analysis} ./logs
for stage in density equil prod; do
  if (( SYMLINK_BULK )); then
    mkdir -p "${SCRATCH_DIR}/${stage}"
    ln -s "${SCRATCH_DIR}/${stage}" "${OUTDIR}/${stage}"
  else
    mkdir -p "${OUTDIR}/${stage}"
  fi
done

[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }

# Record-keeping: preserve the exact input structure alongside the outputs.
echo "[CMD] cp $PDB_IN ${OUTDIR}/input_structure.pdb"
cp "$PDB_IN" "${OUTDIR}/input_structure.pdb"

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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
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
  -ntomp "$NTOMP_SERIAL" \
  > "${EM_DIR}/mdrun_em.log" 2>&1

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
    GEN_VEL="yes"
    GEN_TEMP="${T_MIN}"
    CONT_FLAGS="gen-vel = yes
gen-temp = ${T_MIN}
gen-seed = ${SEED}"
    SEG_DT="0.001"
    SEG_STEPS="${DENSITY_SEG_STEPS}"
  else
    CONT_FLAGS="gen-vel = no
continuation = yes"
    SEG_DT="${DT_PS}"
    SEG_STEPS="${DENSITY_SEG_STEPS}"
  fi

  cat > "${DENSITY_DIR}/density_seg${seg}.mdp" <<EOF
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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
DispCorr    = ${DISPCORR}
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

ld-seed     = ${SEED}

${CONT_FLAGS}

define      = -DPOSRES
EOF

  $GMX grompp \
    -f "${DENSITY_DIR}/density_seg${seg}.mdp" \
    -c "$LAST_GRO" \
    -r "$LAST_GRO" \
    -p "$SYS_TOP" \
    -o "${DENSITY_DIR}/density_seg${seg}.tpr" \
    -po "${DENSITY_DIR}/mdout_density${seg}.mdp" \
    -maxwarn 5 \
    > "${DENSITY_DIR}/grompp_density${seg}.log" 2>&1

  $GMX mdrun \
    -deffnm "${DENSITY_DIR}/density_seg${seg}" \
    -x "${DENSITY_DIR}/${OUTBASE}_density_seg${seg}.xtc" \
    -ntomp "$NTOMP_SERIAL" -nb gpu -pme gpu -bonded gpu \
    > "${DENSITY_DIR}/mdrun_density${seg}.log" 2>&1

  # Extract volume
  echo Volume | $GMX energy \
    -f "${DENSITY_DIR}/density_seg${seg}.edr" \
    -o "${DENSITY_DIR}/volume_seg${seg}.xvg" \
    > /dev/null 2>&1

  vol_avg=$(awk '!/^[#@]/ {sum+=$2; n++} END {if(n>0) printf "%.6f", sum/n}' \
            "${DENSITY_DIR}/volume_seg${seg}.xvg")
  echo "[INFO] Segment $seg average volume: ${vol_avg} nm^3"

  if [[ -n "$prev_vol" ]]; then
    rel_change=$(python3 -c "
prev=$prev_vol; curr=$vol_avg
print(abs(curr-prev)/prev if prev!=0 else 1.0)
")
    echo "[INFO] Relative volume change: $rel_change"

    if (( seg >= DENSITY_MIN_SEG )); then
      converged=$(python3 -c "print(0 if float('$rel_change')>float('${DENSITY_TOL_REL}') else 1)")
      if (( converged )); then
        echo "[INFO] Density converged at segment $seg"
        break
      fi
    fi
  fi

  prev_vol="$vol_avg"
  LAST_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
done

if (( seg == DENSITY_MAX_SEG )); then
  echo "[WARNING] Reached DENSITY_MAX_SEG=${DENSITY_MAX_SEG} without convergence; using last segment"
fi

DENSITY_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
DENSITY_CPT="${DENSITY_DIR}/density_seg${seg}.cpt"
[[ -f "$DENSITY_GRO" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 6 — Prepare per-replica equilibration directories
###############################################################################
CURRENT_STEP=6
echo "===== STEP 6: Preparing per-replica equilibration (${ENSEMBLE}) ====="

EQUIL_DIR="${OUTDIR}/equil"

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  T=${TARR[$i]}
  RSEED=$(rep_seed "$i")
  REPDIR="${EQUIL_DIR}/rep${idx}"
  mkdir -p "$REPDIR"

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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
DispCorr    = ${DISPCORR}
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T} ${T}

${PCOUPL_BLOCK}

constraints = h-bonds
constraint-algorithm = LINCS

ld-seed     = ${RSEED}
gen-vel     = yes
gen-temp    = ${T}
gen-seed    = ${RSEED}

nstcomm     = 500
EOF

  $GMX grompp \
    -f "${REPDIR}/equil.mdp" \
    -c "$DENSITY_GRO" \
    -r "$DENSITY_GRO" \
    -p "$SYS_TOP" \
    -o "${REPDIR}/equil.tpr" \
    -po "${REPDIR}/mdout_equil.mdp" \
    -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  # No per-file symlink needed: equil/ is a folder symlink into scratch, so
  # `-deffnm equil` inside equil/rep*/ already writes equil.xtc straight to scratch.
done
echo "[OK] equilibration .tpr files prepared for all replicas"

###############################################################################
# STEP 7 — Run per-replica equilibration
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Running per-replica equilibration ====="

# Build the -multidir list
MULTIDIR_EQUIL=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  MULTIDIR_EQUIL+="${EQUIL_DIR}/rep${idx} "
done

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir $MULTIDIR_EQUIL \
  -deffnm equil \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_equil.log" 2>&1

# Verify all replicas finished
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${EQUIL_DIR}/rep${idx}/equil.gro" ]] || { echo "[ERROR] equilibration failed for replica $idx"; exit 1; }
done
echo "[OK] per-replica equilibration done"

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
  RSEED=$(rep_seed "$i")
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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
DispCorr    = ${DISPCORR}
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T} ${T}

${PCOUPL_BLOCK}

constraints = h-bonds
constraint-algorithm = LINCS

ld-seed     = ${RSEED}
gen-vel     = no
continuation = yes

nstcomm     = 500
EOF

  $GMX grompp \
    -f "${REPDIR}/remd.mdp" \
    -c "${EQUIL_DIR}/rep${idx}/equil.gro" \
    -t "${EQUIL_DIR}/rep${idx}/equil.cpt" \
    -p "$SYS_TOP" \
    -o "${REPDIR}/remd.tpr" \
    -po "${REPDIR}/mdout_remd.mdp" \
    -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  # No per-file symlink needed: prod/ is a folder symlink into scratch, so
  # `-deffnm remd` inside prod/rep*/ already writes remd.xtc straight to scratch.
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
  -reseed "$SEED" \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_remd.log" 2>&1

# Verify every replica produced a final structure. A nonzero mpirun already aborts
# via set -e + the ERR trap (which preserves scratch, since step 9 >= PRESERVE_FROM_STEP);
# this additionally catches a silent partial failure where mdrun returns 0 but a
# replica wrote no output. Mirrors the per-replica equilibration check in STEP 7.
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${PROD_DIR}/rep${idx}/remd.gro" ]] || { echo "[ERROR] REMD production failed for replica $idx"; exit 1; }
done
echo "[OK] REMD production done"

###############################################################################
# STEP 10 — Export final PDB
###############################################################################
# Trajectories need no finalizing: under the folder-symlink model the .xtc files
# are already at prod/rep*/remd.xtc (real on scratch via the prod/ symlink, or real
# in OUTDIR when SYMLINK_BULK=0). The scratch archive is completed in STEP 13.
CURRENT_STEP=10
echo "===== STEP 10: Export final PDB ====="

# Export final PDB from replica 000. -pbc mol -center -ur compact (centre group
# Protein, output System) — matching the analysis path (fix_PBC.sh). Plain
# -pbc whole only unbreaks molecules and leaves the protein uncentred in the box.
printf "Protein\nSystem\n" | $GMX trjconv \
  -s "${PROD_DIR}/rep000/remd.tpr" \
  -f "${PROD_DIR}/rep000/remd.gro" \
  -o "${OUTDIR}/${OUTBASE}_final_rep000.pdb" \
  -pbc mol -center -ur compact \
  > /dev/null 2>&1 || true

echo "[OK] Final PDB exported"

###############################################################################
# STEP 11 — Write parameters log
###############################################################################
CURRENT_STEP=11
ELAPSED_S=$(( $(date +%s) - CURRENT_TIME ))

# Write to OUTDIR always; also to scratch when offloading, so each scratch archive
# is self-describing. (SYMLINK_BULK=0 → SCRATCH_DIR empty → OUTDIR-only.)
# ── Provenance capture ───────────────────────────────────────────────────────
# What code and what toolchain actually ran. Resolved HERE so a finished job is
# self-describing on its own — reconstructing a run must not depend on still having
# the submit script, or on the repo being unchanged since.
# Probed once (the heredoc below runs per destination) and fail-soft: a version
# probe must never take the job down, so "did it run?" is kept separate from "what
# did it say?" (GOTCHAS: `VAR=$(cmd | grep ...)` under `set -e`).
# The engine path is spelled out rather than taken from BASH_SOURCE — sbatch runs
# from a temp copy, so BASH_SOURCE does not point at the repo (GOTCHAS).
set +e
GMX_VERSION="$($GMX --version 2>/dev/null | sed -n 's/^GROMACS version:[[:space:]]*//p' | head -1)"
SITE_CONFIG_PATH="$(cd "${GROMACS_SCRIPTS_DIR}/../.." && pwd)/site_config.sh"
set -e
GMX_VERSION="${GMX_VERSION:-unknown}"

PARAM_DESTS=( "${OUTDIR}/parameters.txt" )
[[ -n "$SCRATCH_DIR" ]] && PARAM_DESTS+=( "${SCRATCH_DIR}/parameters.txt" )
for dest in "${PARAM_DESTS[@]}"; do
  cat > "$dest" <<EOF
Elapsed seconds: $ELAPSED_S
Job Name: ${SLURM_JOB_NAME:-remd}
Job ID: ${SLURM_JOB_ID:-N/A}
Replicas: $REPLICAS
Total ns per replica: $TOTAL_NS
Temperature range: ${T_MIN}–${T_MAX} K
Force field: $FF${FF_ALIAS_NOTE}
Water model: $WATER
Box shape: $BOX_SHAPE
Box buffer (nm): $BOX_BUFFER
Neutralize: $NEUTRALIZE
Salt molarity: $SALT_MOLAR
Cutoff (nm): $CUTOFF_NM
tau_t (ps): $TAU_T
Time step (ps): $DT_PS
RNG seed: $SEED$( [[ "$SEED" == "-1" ]] && echo " (random per run)" || echo " (fixed; per-replica = SEED+i)" )
Serial-step OMP threads (EM/density): $NTOMP_SERIAL
Exchange interval (ps): $REPLEX_PS
Exchange interval (steps): $REPLEX_STEPS
Production steps: $PROD_NSTEPS
Production ensemble: $ENSEMBLE$( [[ "$ENSEMBLE" == NPT ]] && echo " (C-rescale, ref-p ${REF_P} bar, tau-p ${TAU_P} ps)" )
Per-replica equil ns: $EQUIL_NS
Output directory: $OUTDIR
Scratch directory: ${SCRATCH_DIR:-(none — SYMLINK_BULK=0, bulk real in OUTDIR)}
GPUs: ${SLURM_GPUS_ON_NODE:-unknown}
Density convergence segments: $seg

# ── Provenance ──
Engine: ${GROMACS_SCRIPTS_DIR}/REMD-gromacs.sbatch
Site config: $SITE_CONFIG_PATH
GROMACS: $GMX_VERSION
GMXRC: $GMXRC
GMXLIB: ${GMXLIB:-(unset)}
Input structure: $PDB_IN
Submitted from: ${SLURM_SUBMIT_DIR:-N/A}
EOF
done

###############################################################################
# STEP 12 — Post-analysis (acceptance rates + PBC/strip/align + RMSD/Rg/RMSF/DSSP
#           + conformational clustering [+ per-chain/inter-chain for complexes], rep000)
###############################################################################
CURRENT_STEP=12
echo "===== STEP 12: Post-analysis (rep000 / 300 K) ====="

ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"

# The whole post-analysis lives in run_analysis.sh (auto-detects MD vs REMD from
# the job layout — here it finds prod/rep000/remd.tpr and also runs the REMD-only
# acceptance-rate QC — activates the Python env, and echoes each [CMD]) so it can
# be re-run by hand later without resubmitting the simulation:
#   bash scripts/analysis/run_analysis.sh "$OUTDIR" [REP]
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR 000"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" 000 \
  || echo "[WARN] run_analysis.sh failed — re-run the command above"

###############################################################################
# STEP 13 — Make the scratch archive self-contained (SYMLINK_BULK=1 only)
###############################################################################
# The bulk stages are already real on scratch (written there through the folder
# symlinks). Copy the small real OUTDIR dirs (analysis/ logs/ em/ build/ +
# parameters.txt + final PDB) into scratch so it stands alone. No-op when
# SYMLINK_BULK=0 (nothing was offloaded).
CURRENT_STEP=13
if [[ -n "$SCRATCH_DIR" ]]; then
  echo "===== STEP 13: Finalizing scratch archive ====="
  echo "[CMD] bash ${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh $OUTDIR $SCRATCH_DIR"
  bash "${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh" "$OUTDIR" "$SCRATCH_DIR" \
    || echo "[WARN] finalize_outputs.sh failed — scratch archive may be incomplete; re-run the command above"
else
  echo "===== STEP 13: SYMLINK_BULK=0 — bulk is real in OUTDIR; no scratch archive to finalize ====="
fi

###############################################################################
# STEP 14 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:   $OUTDIR"
if [[ -n "$SCRATCH_DIR" ]]; then
  echo "[INFO] Stage dirs:         prod/ equil/ density/ are symlinks into scratch"
  echo "[INFO] Scratch archive:    $SCRATCH_DIR  (self-contained: bulk + copied analysis/logs/em/build)"
else
  echo "[INFO] Stage dirs:         prod/ equil/ density/ are real in OUTDIR (SYMLINK_BULK=0)"
fi
echo "[INFO] REMD exchange log:  ${PROD_DIR}/rep000/remd.log  (or check all rep dirs)"
echo "[INFO] Acceptance rates:   ${OUTDIR}/analysis/remd_acceptance.csv"
echo "[INFO] Processed traj:     ${OUTDIR}/analysis/remd_rep000_stripped_aligned.xtc"
echo "=========================================="
