#!/bin/bash
#SBATCH --job-name=rest2
#SBATCH -p pi_keating
#SBATCH -C rocky8
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --gres=gpu:l40s:2
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
# The plumed CLI pulls in libhcoll/libocoms transitively (openmpi runtime); some GPU nodes
# lack /opt/mellanox, so prepend the staged compat copies. (This 2023.5 gmx_mpi build does
# NOT link them — unlike the 2024.3 build — but the prepend is a harmless no-op for it.)
export LD_LIBRARY_PATH="${HCOLL_COMPAT_DIR}:${LD_LIBRARY_PATH:-}"
GMX="${GMX:-gmx_mpi}"

command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }

# ---- CUDA MPS: let the multiple ranks sharing each GPU run kernels concurrently. ----
# A/B benchmarked +13.7% aggregate throughput on 16-rep REST2 / 4x L40S (run-to-run
# noise 0.5%) -> ~12% wall-clock per production run, for free. This engine is nodes=1,
# so ONE daemon on the launch node covers every rank. If ever made multi-node, MPS must
# be started per-node via srun instead. Teardown is via the EXIT trap below (fires on
# normal completion AND on the ERR path, since cleanup_scratch_on_error ends in exit);
# a hard SIGKILL (e.g. SLURM timeout) leaves the daemon to the job cgroup to reap.
start_mps() {
  command -v nvidia-cuda-mps-control >/dev/null 2>&1 \
    || { echo "[WARN] nvidia-cuda-mps-control not found — running WITHOUT MPS"; return 0; }
  export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-${SLURM_JOB_ID:-$$}"
  export CUDA_MPS_LOG_DIRECTORY="/tmp/nvidia-mps-log-${SLURM_JOB_ID:-$$}"
  mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
  nvidia-cuda-mps-control -d \
    && echo "[OK] CUDA MPS daemon started ($CUDA_MPS_PIPE_DIRECTORY)" \
    || { echo "[WARN] MPS daemon failed to start — continuing without it"; \
         unset CUDA_MPS_PIPE_DIRECTORY CUDA_MPS_LOG_DIRECTORY; }
}
stop_mps() {
  [[ -n "${CUDA_MPS_PIPE_DIRECTORY:-}" ]] || return 0
  echo quit | nvidia-cuda-mps-control 2>/dev/null || true
  echo "[OK] CUDA MPS daemon stopped"
  unset CUDA_MPS_PIPE_DIRECTORY CUDA_MPS_LOG_DIRECTORY
  return 0
}
trap stop_mps EXIT
start_mps
# Probe the binary for -hrex support. Run mdrun -h ONCE, capture its combined output, and
# check that gmx_mpi actually RAN (rc==0) BEFORE interpreting the result. If we instead piped
# straight into `grep` (as before), a binary that fails to launch — a bad MPI/CUDA/library env
# on a compute node — produces zero 'hrex' matches, grep exits 1, and under `set -e` the
# assignment silently kills the whole job with an empty log (only "JOB START" printed). The
# `2>&1 | grep` also swallowed gmx's real error message. So: disable set -e around the probe,
# inspect the exit code explicitly, and surface the actual failure. Fail loud, never silent.
set +e
_MDRUN_HELP="$("$GMX" mdrun -h 2>&1)"
_MDRUN_RC=$?
set -e
if (( _MDRUN_RC != 0 )); then
  echo "[ERROR] '$GMX mdrun -h' exited $_MDRUN_RC — the GROMACS binary could not run on this node."
  echo "        This is an environment problem (MPI / CUDA / shared libraries), NOT a missing"
  echo "        hrex feature. Real output from the binary follows:"
  echo "-------------------------------------------------------------------------------"
  echo "$_MDRUN_HELP"
  echo "-------------------------------------------------------------------------------"
  exit 1
fi
grep -qi 'hrex' <<<"$_MDRUN_HELP" \
  || { echo "[ERROR] this GROMACS has no -hrex (not PLUMED-patched for hrex)"; exit 1; }
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
# Output model: 1 (default) = bulk stage dirs (density/ equil/ prod/) are folder
# symlinks into scratch (quota-safe, laptop-sync-friendly); 0 = everything real in
# OUTDIR with no scratch offload. (topol/ stays real either way — small + a run record.)
SYMLINK_BULK="${SYMLINK_BULK:-1}"

FF="${FF:-amber99sb-ildn}"
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

# REST2 scales the SOLUTE Hamiltonian via `plumed partial_tempering`, which rewrites
# [ atoms ] [ atomtypes ] [ nonbond_params ] [ pairtypes ] [ dihedrals ] — and has NO
# handling for [ cmap ] / [ cmaptypes ]. Measured: those sections are byte-identical
# between lambda=1.0 and lambda=0.5. CHARMM's CMAP is a protein-backbone cross-term, so
# it sits entirely inside the hot region: a CHARMM REST2 run would scale charges/LJ/
# dihedrals but leave the backbone conformational term at full strength — a partially
# scaled Hamiltonian that neither the scale=1.0 sanity pair nor the hrex acceptance gate
# can detect. Fail here rather than produce a plausible-looking wrong ensemble.
# AMBER (amber99sb-ildn, amber14sb) has no CMAP and is unaffected. T-REMD and plain MD
# scale nothing, so CHARMM is fully correct there. See knowledgebase/GOTCHAS.md.
if [[ "$FF" == charmm* ]]; then
  echo "[ERROR] FF='${FF}' is not supported by the REST2 engine." >&2
  echo "[ERROR] plumed partial_tempering does not scale CHARMM CMAP ([ cmap ]/[ cmaptypes ])," >&2
  echo "[ERROR] so the solute Hamiltonian would be only PARTIALLY scaled and the run would" >&2
  echo "[ERROR] look healthy while sampling the wrong ensemble. Use T-REMD (REMD-gromacs.sbatch)" >&2
  echo "[ERROR] for CHARMM enhanced sampling, or run REST2 with an AMBER force field." >&2
  exit 1
fi
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
GAMMA_LN="${GAMMA_LN:-2.0}"

# RNG seed. -1 (default) = random per run; SEED>=0 pins velocity generation, the stochastic
# V-rescale/C-rescale streams, and the replica-exchange RNG. Per-replica stages offset it by
# the replica index (SEED+i). NOTE: GPU runs are not bit-for-bit reproducible regardless.
SEED="${SEED:--1}"
rep_seed() { if [[ "$SEED" == "-1" ]]; then echo -1; else echo $(( SEED + $1 )); fi; }

REPLEX_PS="${REPLEX_PS:-1.0}"            # exchange attempt interval (ps); keep >= 1 (GPU deadlock gotcha)
TOTAL_NS="${TOTAL_NS:-20}"

HEAT_NS="${HEAT_NS:-0.2}"                # NVT thermalization length (ns); velocities are generated here
EQUIL_NS="${EQUIL_NS:-0.2}"
DENSITY_SEG_STEPS="${DENSITY_SEG_STEPS:-10000}"
DENSITY_MIN_SEG="${DENSITY_MIN_SEG:-8}"
DENSITY_MAX_SEG="${DENSITY_MAX_SEG:-20}"
DENSITY_TOL_REL="${DENSITY_TOL_REL:-0.005}"

ENSEMBLE="${ENSEMBLE:-NVT}"
# Barostat settings — see the REMD engine for the full note. In short: the density stage
# is pressure-coupled under BOTH ensembles, so REF_P is the reference pressure of the box
# production inherits; TAU_P=2.0 keeps C-rescale off the instantaneous virial noise.
REF_P="${REF_P:-1.0}"                    # bar
TAU_P="${TAU_P:-2.0}"                    # ps

export PDB_IN OUTBASE OUTDIR FF WATER BOX_SHAPE BOX_BUFFER NEUTRALIZE \
       SALT_MOLAR REPLICAS NTOMP_SERIAL DT_PS T_MIN T_MAX TEMPS_LIST CUTOFF_NM GAMMA_LN \
       REPLEX_PS TOTAL_NS HEAT_NS EQUIL_NS DENSITY_SEG_STEPS DENSITY_MIN_SEG DENSITY_MAX_SEG \
       DENSITY_TOL_REL ENSEMBLE REF_P TAU_P PRESERVE_SCRATCH_FROM SYMLINK_BULK SEED CFG_REPLICAS
python3 "$VALIDATE_PARAMS" --engine rest2

TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")
REPLEX_STEPS=$(python3 -c "print(int(round(float('${REPLEX_PS}')/float('${DT_PS}'))))")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")
EQUIL_NSTEPS=$(python3 -c "print(int(round(float('${EQUIL_NS}')*1000/float('${DT_PS}'))))")
HEAT_NSTEPS=$(python3 -c "print(int(round(float('${HEAT_NS}')*1000/float('${DT_PS}'))))")
PREFLIGHT_NSTEPS=$(( REPLEX_STEPS * 20 ))   # ~20 exchange attempts for the acceptance gate

case "$PRESERVE_SCRATCH_FROM" in
  always)  PRESERVE_FROM_STEP=0 ;;
  density) PRESERVE_FROM_STEP=6 ;;
  prod)    PRESERVE_FROM_STEP=12 ;;
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
    # PRESERVE_SCRATCH_FROM decides the fate of the scratch dir, full stop — an explicit
    # policy is not overridden here, and "never" means never. Constraint-failure dumps
    # (step*.pdb, in the cwd = SCRATCH_DIR since EM) are only 75 MB each next to a scratch
    # tree that can be many GB, and on a strained structure EM dumps on EVERY run, so
    # keeping scratch whenever dumps exist would quietly disable the knob for that whole
    # class of input. Instead: always write the few-KB report into OUTDIR so the failure
    # stays on record, and let the policy govern the bulky PDBs.
    if (( CURRENT_STEP >= PRESERVE_FROM_STEP )); then
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); PRESERVING scratch: $SCRATCH_DIR"
      bash "${GROMACS_SCRIPTS_DIR}/report_crash_dumps.sh" "$SCRATCH_DIR" "$OUTDIR" retained || true
    else
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); cleaning scratch: $SCRATCH_DIR"
      bash "${GROMACS_SCRIPTS_DIR}/report_crash_dumps.sh" "$SCRATCH_DIR" "$OUTDIR" deleted || true
      rm -rf "$SCRATCH_DIR"
    fi
  fi
  exit $rc
}
trap 'cleanup_scratch_on_error' ERR

# Check disk space where the bulk will actually land (scratch, or OUTDIR when SYMLINK_BULK=0)
SCRATCH_AVAIL=$(df -BG "$DF_TARGET" | awk 'NR==2 {print $4}' | sed 's/G//')
if (( SCRATCH_AVAIL < 100 )); then
  echo "[ERROR] Low disk space in ${DF_TARGET}: ${SCRATCH_AVAIL}GB available (need ~100GB)"; exit 1
fi

# ── Output directory tree (folder-symlink model) ─────────────────────────────
# Small, laptop-worthy dirs stay REAL in OUTDIR (build/ em/ topol/ logs/ analysis/).
# The bulk stage dirs (density/ equil/ prod/) are SYMLINKS into scratch when
# SYMLINK_BULK=1 — created before the run so mdrun writes straight onto scratch and
# the tight-quota pool is untouched — or REAL dirs in OUTDIR when SYMLINK_BULK=0.
# Either way the stage-dir variables point at ${OUTDIR}/<stage>, so nothing
# downstream changes. build/ stays real so the solvate/genion topology rename()
# stays same-filesystem (see gotcha).
mkdir -p "${OUTDIR}"/{build,em,topol,logs,analysis} ./logs
for stage in heat density equil prod; do
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

# ── Where GROMACS drops its crash dumps ──────────────────────────────────────
# On a constraint failure mdrun writes step<N>b.pdb / step<N>c.pdb (coordinates before
# and after constraining) into the CURRENT WORKING DIRECTORY — which for an sbatch is the
# submit directory, not the job's OUTDIR. On a large system those are ~75 MB *each*, so a
# strained input silently drops hundreds of MB next to the submit scripts on the
# tight-quota pool. Move cwd to where the bulk already goes, so dumps land with the job.
#
# THIS MUST COME AFTER THE BUILD STEP, NOT BEFORE IT. solvate/genion finish their topology
# edit with a rename() that is fatal across filesystems, so while they run the cwd has to
# stay on the same filesystem as the topology (see the OUTDIR-filesystem gotcha). Moving
# this above the build makes STEP 3 die in solvate. Every mdrun runs at or after EM, so
# here is early enough to catch them all.
#
# Otherwise safe: every path used below is absolute; SLURM opened its -o/-e files relative
# to the submit dir before this script body started; $SLURM_SUBMIT_DIR is an env var; and
# multi-replica mdrun (-multidir) chdirs each rank into its own rep dir regardless.
cd "${SCRATCH_DIR:-$OUTDIR}"

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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
pbc         = xyz
define      = -DPOSRES
EOF
$GMX grompp -f "${EM_DIR}/em.mdp" -c "$SYS_GRO" -r "$SYS_GRO" -p "$SYS_TOP" \
  -o "${EM_DIR}/em.tpr" -po "${EM_DIR}/mdout_em.mdp" -maxwarn 5 > "${EM_DIR}/grompp_em.log" 2>&1
$GMX mdrun -v -deffnm "${EM_DIR}/em" -ntomp "$NTOMP_SERIAL" > "${EM_DIR}/mdrun_em.log" 2>&1
[[ -f "${EM_DIR}/em.gro" ]] || { echo "[ERROR] Minimization failed"; exit 1; }
echo "[OK] Minimization complete"

###############################################################################
# STEP 5 — Heat (NVT, restrained, thermalize at T_MIN)
###############################################################################
# Thermalize at constant volume BEFORE the barostat is switched on, matching the MD and
# T-REMD engines. Previously the density stage generated velocities itself and ran its
# first segment at dt=0.001 as a crash guard; that put C-rescale on a minimized
# configuration whose instantaneous virial pressure is far from equilibrium. Runs on the
# UNSCALED topology — solute tempering starts at STEP 7.
CURRENT_STEP=5
echo "===== STEP 5: Heat (NVT, restrained) ====="
HEAT_DIR="${OUTDIR}/heat"
cat > "${HEAT_DIR}/heat.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${HEAT_NSTEPS}
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
ref-t       = ${T_MIN} ${T_MIN}
pcoupl      = no
constraints = h-bonds
constraint-algorithm = LINCS
ld-seed     = ${SEED}
gen-vel     = yes
gen-temp    = ${T_MIN}
gen-seed    = ${SEED}
nstcomm     = 500
define      = -DPOSRES
EOF
echo "[CMD] $GMX grompp -f ${HEAT_DIR}/heat.mdp -c ${EM_DIR}/em.gro -r ${EM_DIR}/em.gro -p $SYS_TOP -o ${HEAT_DIR}/heat.tpr"
$GMX grompp -f "${HEAT_DIR}/heat.mdp" -c "${EM_DIR}/em.gro" -r "${EM_DIR}/em.gro" -p "$SYS_TOP" \
  -o "${HEAT_DIR}/heat.tpr" -po "${HEAT_DIR}/mdout_heat.mdp" -maxwarn 5 \
  > "${HEAT_DIR}/grompp.log" 2>&1
echo "[CMD] $GMX mdrun -deffnm ${HEAT_DIR}/heat -ntomp $NTOMP_SERIAL -nb gpu -pme gpu -bonded gpu"
$GMX mdrun -deffnm "${HEAT_DIR}/heat" -x "${HEAT_DIR}/heat.xtc" \
  -ntomp "$NTOMP_SERIAL" -nb gpu -pme gpu -bonded gpu > "${OUTDIR}/logs/mdrun_heat.log" 2>&1
[[ -f "${HEAT_DIR}/heat.gro" && -f "${HEAT_DIR}/heat.cpt" ]] || { echo "[ERROR] Heat step failed (no heat.gro/heat.cpt). The density stage resumes from heat.cpt, so HEAT_NS must be > 0."; exit 1; }
echo "[OK] Heat step done"

###############################################################################
# STEP 6 — NPT density equilibration (iterative convergence)
###############################################################################
CURRENT_STEP=6
echo "===== STEP 6: NPT density equilibration (iterative) ====="
DENSITY_DIR="${OUTDIR}/density"
# Fixed position-restraint reference: grompp -r sets the coordinates the harmonic
# restraints pull toward, so passing the previous segment's output (as this loop used to)
# re-anchors the restraint to wherever the protein has already drifted — resisting
# per-segment displacement but never CUMULATIVE displacement.
POSRES_REF="${EM_DIR}/em.gro"
LAST_GRO="${HEAT_DIR}/heat.gro"
LAST_CPT="${HEAT_DIR}/heat.cpt"
volumes=()
seg=0
while (( seg < DENSITY_MAX_SEG )); do
  seg=$((seg+1))
  echo "[INFO] NPT segment $seg"
  cat > "${DENSITY_DIR}/density_seg${seg}.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${DENSITY_SEG_STEPS}
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
tau-p       = ${TAU_P}
ref-p       = ${REF_P}
compressibility = 4.5e-5
refcoord-scaling = com
constraints = h-bonds
constraint-algorithm = LINCS
ld-seed     = ${SEED}
gen-vel     = no
continuation = yes
define      = -DPOSRES
EOF
  # -t carries the integrator state (V-rescale + C-rescale internal state and their RNG
  # streams) across the segment boundary; a .gro holds only positions, velocities and box,
  # so without it the barostat restarts cold each segment and kicks the volume — the very
  # quantity the convergence test measures.
  $GMX grompp -f "${DENSITY_DIR}/density_seg${seg}.mdp" -c "$LAST_GRO" -r "$POSRES_REF" -t "$LAST_CPT" -p "$SYS_TOP" \
    -o "${DENSITY_DIR}/density_seg${seg}.tpr" -po "${DENSITY_DIR}/mdout_density${seg}.mdp" -maxwarn 5 \
    > "${DENSITY_DIR}/grompp_density${seg}.log" 2>&1
  $GMX mdrun -deffnm "${DENSITY_DIR}/density_seg${seg}" \
    -x "${DENSITY_DIR}/${OUTBASE}_density_seg${seg}.xtc" \
    -ntomp "$NTOMP_SERIAL" -nb gpu -pme gpu -bonded gpu > "${DENSITY_DIR}/mdrun_density${seg}.log" 2>&1
  echo Volume | $GMX energy -f "${DENSITY_DIR}/density_seg${seg}.edr" \
    -o "${DENSITY_DIR}/volume_seg${seg}.xvg" > /dev/null 2>&1
  vol_avg=$(awk '!/^[#@]/ {sum+=$2; n++} END {if(n>0) printf "%.6f", sum/n}' "${DENSITY_DIR}/volume_seg${seg}.xvg")
  # A missing/empty volume would otherwise reach density_converged.py as an empty argv
  # entry and die there as an argparse error about a float, several layers from the cause.
  [[ -n "$vol_avg" ]] || { echo "[ERROR] no Volume in ${DENSITY_DIR}/volume_seg${seg}.xvg — did 'gmx energy' fail for segment ${seg}?"; exit 1; }
  echo "[INFO] Segment $seg average volume: ${vol_avg} nm^3"
  volumes+=("$vol_avg")
  # Plateau test over the trailing DENSITY_MIN_SEG segments, not a consecutive-segment
  # difference — a steady drift just under DENSITY_TOL_REL per segment passes the latter on
  # every comparison while the box keeps contracting. See density_converged.py.
  conv_line=$(python3 "${GROMACS_SCRIPTS_DIR}/density_converged.py" \
                --tol "$DENSITY_TOL_REL" --min-seg "$DENSITY_MIN_SEG" "${volumes[@]}")
  converged="${conv_line%% *}"            # leading verdict: 1 = converged, 0 = not
  echo "[INFO] density: ${conv_line#* }"  # trailing reason, next to the volumes above
  if (( converged )); then echo "[INFO] Density converged at segment $seg"; break; fi
  LAST_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
  LAST_CPT="${DENSITY_DIR}/density_seg${seg}.cpt"
done
if (( ! converged )); then
  echo "[WARNING] Reached DENSITY_MAX_SEG=${DENSITY_MAX_SEG} without convergence; using last segment"
fi
DENSITY_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
DENSITY_CPT="${DENSITY_DIR}/density_seg${seg}.cpt"
[[ -f "$DENSITY_GRO" && -f "$DENSITY_CPT" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 7 — Scaled topologies (grompp -pp → mark protein → partial_tempering)
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Generating per-replica scaled topologies ====="
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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
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
# STEP 8 — Prepare per-replica equilibration (all at T_MIN, scaled Hamiltonians)
###############################################################################
CURRENT_STEP=8
echo "===== STEP 8: Preparing per-replica equilibration (${ENSEMBLE}) ====="
EQUIL_DIR="${OUTDIR}/equil"
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  RSEED=$(rep_seed "$i")
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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
DispCorr    = ${DISPCORR}
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_MIN} ${T_MIN}
${PCOUPL_BLOCK}
constraints = h-bonds
constraint-algorithm = LINCS
ld-seed     = ${RSEED}
gen-vel     = yes
gen-temp    = ${T_MIN}
gen-seed    = ${RSEED}
nstcomm     = 500
EOF
  $GMX grompp -f "${REPDIR}/equil.mdp" -c "$DENSITY_GRO" -r "$DENSITY_GRO" \
    -p "${TOPOL_DIR}/topol_rep${idx}.top" -o "${REPDIR}/equil.tpr" -po "${REPDIR}/mdout_equil.mdp" -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  # No per-file symlink needed: equil/ is a folder symlink into scratch (or real when
  # SYMLINK_BULK=0), so `-deffnm equil` inside equil/rep*/ writes equil.xtc there directly.
done
echo "[OK] equilibration .tpr files prepared"

###############################################################################
# STEP 9 — Run per-replica equilibration (independent MD, NO exchange)
###############################################################################
CURRENT_STEP=9
echo "===== STEP 9: Running per-replica equilibration ====="
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
# STEP 10 — Prepare REST2 production inputs
###############################################################################
CURRENT_STEP=10
echo "===== STEP 10: Preparing REST2 production inputs ====="
PROD_DIR="${OUTDIR}/prod"
NSTXOUT_COMP=$(python3 -c "print(int(${REPLEX_STEPS}*10))")
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  RSEED=$(rep_seed "$i")
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
${VDW_BLOCK}rvdw        = ${CUTOFF_NM}
DispCorr    = ${DISPCORR}
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_MIN} ${T_MIN}
${PCOUPL_BLOCK}
constraints = h-bonds
constraint-algorithm = LINCS
ld-seed     = ${RSEED}
gen-vel     = no
continuation = yes
nstcomm     = 500
EOF
  $GMX grompp -f "${REPDIR}/rest2.mdp" -c "${EQUIL_DIR}/rep${idx}/equil.gro" -t "${EQUIL_DIR}/rep${idx}/equil.cpt" \
    -p "${TOPOL_DIR}/topol_rep${idx}.top" -o "${REPDIR}/rest2.tpr" -po "${REPDIR}/mdout_rest2.mdp" -maxwarn 5 \
    > "${REPDIR}/grompp.log" 2>&1
  printf "# REST2 hrex, no CVs (rep%s)\n" "$idx" > "${REPDIR}/plumed.dat"
  # No per-file symlink needed: prod/ is a folder symlink into scratch (or real when
  # SYMLINK_BULK=0), so `-deffnm rest2` inside prod/rep*/ writes rest2.xtc there directly.
done
echo "[OK] REST2 .tpr files prepared"

###############################################################################
# STEP 11 — hrex pre-flight acceptance gate (FAIL LOUD if zero exchanges)
###############################################################################
CURRENT_STEP=11
echo "===== STEP 11: hrex pre-flight acceptance check (~20 exchange attempts) ====="
# Ephemeral workspace (rm -rf'd below): on scratch when offloading, else under OUTDIR.
PRE_DIR="${SCRATCH_DIR:-$OUTDIR}/preflight"
MULTIDIR_PRE=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  mkdir -p "${PRE_DIR}/rep${idx}"
  cp "${PROD_DIR}/rep${idx}/rest2.tpr" "${PROD_DIR}/rep${idx}/plumed.dat" "${PRE_DIR}/rep${idx}/"
  MULTIDIR_PRE+="${PRE_DIR}/rep${idx} "
done
mpirun -np "$REPLICAS" "$GMX" mdrun -multidir $MULTIDIR_PRE -deffnm rest2 \
  -plumed plumed.dat -hrex -replex "$REPLEX_STEPS" -reseed "$SEED" -nsteps "$PREFLIGHT_NSTEPS" \
  -ntomp 1 -nb gpu -pme gpu -bonded gpu -pin on -dlb no > "${OUTDIR}/logs/mdrun_preflight.log" 2>&1
echo "[CMD] python3 ${GROMACS_SCRIPTS_DIR}/check_hrex_acceptance.py ${PRE_DIR}/rep000/rest2.log"
python3 "${GROMACS_SCRIPTS_DIR}/check_hrex_acceptance.py" "${PRE_DIR}/rep000/rest2.log"
rm -rf "$PRE_DIR"
echo "[OK] hrex is exchanging — proceeding to production"

###############################################################################
# STEP 12 — Run REST2 production
###############################################################################
CURRENT_STEP=12
echo "===== STEP 12: Running REST2 production ====="
echo "[INFO] ${REPLICAS} replicas × ${TOTAL_NS} ns, exchange every ${REPLEX_PS} ps, ${T_MIN} K physical"
MULTIDIR_PROD=""
for i in $(seq 0 $((REPLICAS-1))); do MULTIDIR_PROD+="${PROD_DIR}/rep$(printf %03d $i) "; done
mpirun -np "$REPLICAS" "$GMX" mdrun -multidir $MULTIDIR_PROD -deffnm rest2 \
  -plumed plumed.dat -hrex -replex "$REPLEX_STEPS" -reseed "$SEED" \
  -ntomp 1 -nb gpu -pme gpu -bonded gpu -pin on -dlb no > "${OUTDIR}/logs/mdrun_rest2.log" 2>&1
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${PROD_DIR}/rep${idx}/rest2.gro" ]] || { echo "[ERROR] REST2 production failed for replica $idx"; exit 1; }
done
echo "[OK] REST2 production done"

###############################################################################
# STEP 13 — Export final PDB
###############################################################################
# Trajectories need no finalizing: under the folder-symlink model the .xtc files
# are already at prod/rep*/rest2.xtc (real on scratch via the prod/ symlink, or real
# in OUTDIR when SYMLINK_BULK=0). The scratch archive is completed in STEP 16.
CURRENT_STEP=13
echo "===== STEP 13: Export final PDB ====="
# -pbc mol -center -ur compact (centre group Protein, output System) — matching the
# analysis path (fix_PBC.sh). Plain -pbc whole only unbreaks molecules and leaves
# the protein uncentred in the box.
printf "Protein\nSystem\n" | $GMX trjconv -s "${PROD_DIR}/rep000/rest2.tpr" -f "${PROD_DIR}/rep000/rest2.gro" \
  -o "${OUTDIR}/${OUTBASE}_final_rep000.pdb" -pbc mol -center -ur compact > /dev/null 2>&1 || true
echo "[OK] Final PDB exported"

###############################################################################
# STEP 14 — Write parameters log
###############################################################################
CURRENT_STEP=14
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
PLUMED_VERSION="$(plumed info --version 2>/dev/null)"
set -e
GMX_VERSION="${GMX_VERSION:-unknown}"
PLUMED_VERSION="${PLUMED_VERSION:-unknown}"

PARAM_DESTS=( "${OUTDIR}/parameters.txt" )
[[ -n "$SCRATCH_DIR" ]] && PARAM_DESTS+=( "${SCRATCH_DIR}/parameters.txt" )
for dest in "${PARAM_DESTS[@]}"; do
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
Force field: $FF${FF_ALIAS_NOTE}
Water model: $WATER
Box shape: $BOX_SHAPE ; buffer (nm): $BOX_BUFFER
Neutralize: $NEUTRALIZE ; Salt molarity: $SALT_MOLAR
Cutoff (nm): $CUTOFF_NM ; tau_t (ps): $TAU_T ; Time step (ps): $DT_PS
RNG seed: $SEED$( [[ "$SEED" == "-1" ]] && echo " (random per run)" || echo " (fixed; per-replica = SEED+i)" )
Exchange interval: ${REPLEX_PS} ps (${REPLEX_STEPS} steps)
Production ensemble: $ENSEMBLE$( [[ "$ENSEMBLE" == NPT ]] && echo " (C-rescale, ref-p ${REF_P} bar, tau-p ${TAU_P} ps)" )
Heat ns: $HEAT_NS
Per-replica equil ns: $EQUIL_NS
Output directory: $OUTDIR
Scratch directory: ${SCRATCH_DIR:-(none — SYMLINK_BULK=0, bulk real in OUTDIR)}
GPUs: ${SLURM_GPUS_ON_NODE:-unknown}
Density convergence segments: $seg$( (( converged )) && echo " (plateau reached)" || echo " (DENSITY_MAX_SEG cap — NOT converged)" )

# ── Provenance ──
Engine: ${GROMACS_SCRIPTS_DIR}/REST2-gromacs.sbatch
Site config: $SITE_CONFIG_PATH
GROMACS: $GMX_VERSION
REST2_GMXRC: $REST2_GMXRC
PLUMED: $PLUMED_VERSION
GMXLIB: ${GMXLIB:-(unset)}
Input structure: $PDB_IN
Submitted from: ${SLURM_SUBMIT_DIR:-N/A}
EOF
done

###############################################################################
# STEP 15 — Post-analysis (rep000 = lambda=1 = physical T_MIN ensemble)
###############################################################################
CURRENT_STEP=15
echo "===== STEP 15: Post-analysis (rep000 / lambda=1 / ${T_MIN} K) ====="
ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR 000"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" 000 \
  || echo "[WARN] run_analysis.sh failed — re-run the command above"

###############################################################################
# STEP 16 — Make the scratch archive self-contained (SYMLINK_BULK=1 only)
###############################################################################
# The bulk stages are already real on scratch (written there through the folder
# symlinks). Copy the small real OUTDIR dirs (analysis/ logs/ em/ build/ topol/ +
# parameters.txt + final PDB) into scratch so it stands alone. No-op when
# SYMLINK_BULK=0 (nothing was offloaded).
CURRENT_STEP=16
if [[ -n "$SCRATCH_DIR" ]]; then
  echo "===== STEP 16: Finalizing scratch archive ====="
  echo "[CMD] bash ${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh $OUTDIR $SCRATCH_DIR"
  bash "${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh" "$OUTDIR" "$SCRATCH_DIR" \
    || echo "[WARN] finalize_outputs.sh failed — scratch archive may be incomplete; re-run the command above"
else
  echo "===== STEP 16: SYMLINK_BULK=0 — bulk is real in OUTDIR; no scratch archive to finalize ====="
fi

# ── Constraint-failure check ─────────────────────────────────────────────────
# Dumps live on scratch (cwd since EM) so they stay off the tight-quota pool, but
# nobody browses scratch — so surface them in OUTDIR and in this log. Writes
# OUTDIR/CONSTRAINT_FAILURES.txt only when there are dumps, and clears a stale one
# otherwise, so the file's presence always refers to THIS run.
echo "[CMD] bash ${GROMACS_SCRIPTS_DIR}/report_crash_dumps.sh ${SCRATCH_DIR:-$OUTDIR} $OUTDIR"
bash "${GROMACS_SCRIPTS_DIR}/report_crash_dumps.sh" "${SCRATCH_DIR:-$OUTDIR}" "$OUTDIR" || true

###############################################################################
# STEP 17 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:   $OUTDIR"
if [[ -n "$SCRATCH_DIR" ]]; then
  echo "[INFO] Stage dirs:         prod/ equil/ density/ are symlinks into scratch"
  echo "[INFO] Scratch archive:    $SCRATCH_DIR  (self-contained: bulk + copied analysis/logs/em/build/topol)"
else
  echo "[INFO] Stage dirs:         prod/ equil/ density/ are real in OUTDIR (SYMLINK_BULK=0)"
fi
echo "[INFO] rep000 = lambda=1 = ${T_MIN} K physical ensemble (analyze directly)"
echo "[INFO] Exchange log:       ${PROD_DIR}/rep000/rest2.log"
echo "[INFO] Acceptance rates:   ${OUTDIR}/analysis/remd_acceptance.csv"
echo "=========================================="
