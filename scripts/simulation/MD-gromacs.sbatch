#!/bin/bash
#SBATCH --job-name=md
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# GROMACS plain production MD single-node workflow
# ---------------------------------------------------------------------------
# Standard single-temperature MD pipeline (no replica exchange):
#   EM → heat (NVT, restrained) → density (NPT, restrained, iterative)
#      → relax (NPT, unrestrained, optional) → production (NPT) → post-analysis.
# Restraints (-DPOSRES on protein heavy atoms) are held through EM/heat/density
# to preserve the input structure, then released for the unrestrained equilibration
# so production records a clean equilibrium ensemble. Equilibration and production
# share one ensemble (V-rescale + C-rescale).
#
# Intended usage:
#   Copy example/submit_jobs/submit_MD.sh, set your parameters, run it.
#   That script sets GROMACS_SCRIPTS_DIR and passes parameters via --export.
#
# Advanced / direct usage:
#   sbatch --export=ALL,PDB_IN=my.pdb,... MD-gromacs.sbatch
#   sbatch MD-gromacs.sbatch config.sh   # source a config file as $1
#
# NOTE: this build ships only `gmx_mpi` (no `gmx`/`gmx mdrun`). Every step —
# including production — runs as a single-rank `gmx_mpi mdrun` with no mpirun.
###############################################################################

set -euo pipefail

CURRENT_TIME=$(date +%s)
CURRENT_STEP=0
PRESERVE_SCRATCH_FROM="${PRESERVE_SCRATCH_FROM:-prod}"   # prod | density | always | never (resolved to a step index in STEP 1)

# ── Source optional config file ──────────────────────────────────────────────
: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — submit via example/submit_jobs/submit_MD.sh}"
VALIDATE_PARAMS="${GROMACS_SCRIPTS_DIR}/validate_params.py"
if [[ -n "${1:-}" ]]; then
  [[ -f "$1" ]] || { echo "[ERROR] config file not found: $1"; exit 1; }
  python3 "$VALIDATE_PARAMS" --engine md --check-keys "$1"   # reject typo'd keys before sourcing
  echo "[INFO] Sourcing config: $1"
  # shellcheck disable=SC1090
  source "$1"
fi

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
OUTDIR="${OUTDIR:-md_${OUTBASE}}"

# Force field & water model
FF="${FF:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
BOX_SHAPE="${BOX_SHAPE:-dodecahedron}"
BOX_BUFFER="${BOX_BUFFER:-1.0}"          # nm
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

# Simulation parameters
DT_PS="${DT_PS:-0.002}"
T_SIM="${T_SIM:-300}"                     # single production temperature (K)
CUTOFF_NM="${CUTOFF_NM:-1.0}"             # 10 Å = 1.0 nm
GAMMA_LN="${GAMMA_LN:-2.0}"               # ps^-1  (tau_t = 1/gamma)

# Production
TOTAL_NS="${TOTAL_NS:-100}"               # total production (ns)
TRAJ_PS="${TRAJ_PS:-10}"                  # write a trajectory frame every this many ps

# Equilibration
HEAT_NS="${HEAT_NS:-0.2}"                 # NVT heat length (ns)
DENSITY_SEG_STEPS="${DENSITY_SEG_STEPS:-10000}"
DENSITY_MIN_SEG="${DENSITY_MIN_SEG:-8}"
DENSITY_MAX_SEG="${DENSITY_MAX_SEG:-20}"
DENSITY_TOL_REL="${DENSITY_TOL_REL:-0.005}"
RELAX_NS="${RELAX_NS:-0}"               # unrestrained NPT equil after density (ns); 0 = OFF, restraints release at production start

# CPU threads for the single-rank mdrun
NTOMP="${NTOMP:-${SLURM_CPUS_PER_TASK:-8}}"

# Validate resolved parameters before deriving anything or building. Export the job
# params so validate_params.py reads them from the environment (covers config- and
# --export-set values alike). Keep this list in sync with validate_params.py.
export PDB_IN OUTBASE OUTDIR FF WATER BOX_SHAPE BOX_BUFFER NEUTRALIZE SALT_MOLAR \
       DT_PS T_SIM CUTOFF_NM GAMMA_LN TOTAL_NS TRAJ_PS HEAT_NS DENSITY_SEG_STEPS \
       DENSITY_MIN_SEG DENSITY_MAX_SEG DENSITY_TOL_REL RELAX_NS NTOMP PRESERVE_SCRATCH_FROM
python3 "$VALIDATE_PARAMS" --engine md

# Derived step counts (after validation, so a non-numeric input fails with a clear
# message above instead of a raw python traceback here).
TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")
NSTXOUT_COMP=$(python3 -c "print(int(round(float('${TRAJ_PS}')/float('${DT_PS}'))))")
HEAT_NSTEPS=$(python3 -c "print(int(round(float('${HEAT_NS}')*1000/float('${DT_PS}'))))")
RELAX_NSTEPS=$(python3 -c "print(int(round(float('${RELAX_NS}')*1000/float('${DT_PS}'))))")

# Resolve the scratch-preservation stage name to this engine's step index (stage names
# mean the same thing across MD/REMD; the index they map to is engine-specific). The
# validator already restricted PRESERVE_SCRATCH_FROM to {prod,density,always,never}.
case "$PRESERVE_SCRATCH_FROM" in
  always)  PRESERVE_FROM_STEP=0 ;;
  density) PRESERVE_FROM_STEP=5 ;;
  prod)    PRESERVE_FROM_STEP=7 ;;
  never)   PRESERVE_FROM_STEP=999 ;;
  *) echo "[ERROR] invalid PRESERVE_SCRATCH_FROM='$PRESERVE_SCRATCH_FROM'"; exit 1 ;;
esac

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
mkdir -p "${OUTDIR}"/{build,em,heat,density,prod,logs,trajectories} ./logs

[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set  (T=${T_SIM} K, ${TOTAL_NS} ns, ${NTOMP} OMP threads)"

###############################################################################
# STEP 2 — Build system (pdb2gmx → editconf → solvate → genion)
###############################################################################
CURRENT_STEP=2
echo "===== STEP 2: Building system ====="

BUILD="${OUTDIR}/build"

echo "[INFO] Running pdb2gmx (ff=${FF}, water=${WATER})…"
$GMX pdb2gmx \
  -f "$PDB_IN" \
  -o "${BUILD}/${OUTBASE}_processed.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  -i "${BUILD}/${OUTBASE}_posre.itp" \
  -ff "$FF" -water "$WATER" -ignh \
  > "${BUILD}/pdb2gmx.log" 2>&1

echo "[INFO] Setting up ${BOX_SHAPE} box with ${BOX_BUFFER} nm buffer…"
$GMX editconf \
  -f "${BUILD}/${OUTBASE}_processed.gro" \
  -o "${BUILD}/${OUTBASE}_box.gro" \
  -bt "$BOX_SHAPE" -d "$BOX_BUFFER" -c \
  > "${BUILD}/editconf.log" 2>&1

echo "[INFO] Solvating…"
$GMX solvate \
  -cp "${BUILD}/${OUTBASE}_box.gro" \
  -cs spc216.gro \
  -o "${BUILD}/${OUTBASE}_solv.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  > "${BUILD}/solvate.log" 2>&1

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

SYS_GRO="${BUILD}/${OUTBASE}_ions.gro"
SYS_TOP="${BUILD}/${OUTBASE}.top"
[[ -f "$SYS_GRO" && -f "$SYS_TOP" ]] || { echo "[ERROR] System build failed"; exit 1; }
echo "[OK] System built"

###############################################################################
# STEP 3 — Energy minimization
###############################################################################
CURRENT_STEP=3
echo "===== STEP 3: Energy minimization ====="

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
  -ntomp "$NTOMP" \
  > "${EM_DIR}/mdrun_em.log" 2>&1

[[ -f "${EM_DIR}/em.gro" ]] || { echo "[ERROR] Minimization failed"; exit 1; }
echo "[OK] Minimization complete"

###############################################################################
# STEP 4 — Heat (NVT, restrained, heat to T_SIM)
###############################################################################
CURRENT_STEP=4
echo "===== STEP 4: Heat (NVT, restrained) ====="

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
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_SIM} ${T_SIM}

pcoupl      = no

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = yes
gen-temp    = ${T_SIM}

nstcomm     = 500

; Restrain protein heavy atoms while solvent relaxes around the structure
define      = -DPOSRES
EOF

# The heat step continues from the minimized structure (em/em.gro), NOT the un-minimized
# build (SYS_GRO) — feeding SYS_GRO here discards EM and starts the heat step on clashing
# coordinates, which explodes at step 0 (LINCS error / huge +LJ). The position-
# restraint reference (-r) is the same minimized structure.
$GMX grompp \
  -f "${HEAT_DIR}/heat.mdp" \
  -c "${EM_DIR}/em.gro" \
  -r "${EM_DIR}/em.gro" \
  -p "$SYS_TOP" \
  -o "${HEAT_DIR}/heat.tpr" \
  -po "${HEAT_DIR}/mdout_heat.mdp" \
  -maxwarn 5 \
  > "${HEAT_DIR}/grompp.log" 2>&1

$GMX mdrun \
  -deffnm "${HEAT_DIR}/heat" \
  -x "${SCRATCH_DIR}/heat.xtc" \
  -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu -pin on \
  > "${OUTDIR}/logs/mdrun_heat.log" 2>&1

[[ -f "${HEAT_DIR}/heat.gro" ]] || { echo "[ERROR] Heat step failed"; exit 1; }
echo "[OK] Heat step done"

###############################################################################
# STEP 5 — NPT density equilibration (restrained, iterative convergence)
###############################################################################
CURRENT_STEP=5
echo "===== STEP 5: NPT density equilibration (iterative) ====="

DENSITY_DIR="${OUTDIR}/density"
LAST_GRO="${HEAT_DIR}/heat.gro"
LAST_CPT="${HEAT_DIR}/heat.cpt"
prev_vol=""
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
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_SIM} ${T_SIM}

pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 1.0
ref-p       = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = no
continuation = yes

define      = -DPOSRES
EOF

  $GMX grompp \
    -f "${DENSITY_DIR}/density_seg${seg}.mdp" \
    -c "$LAST_GRO" \
    -r "$LAST_GRO" \
    -t "$LAST_CPT" \
    -p "$SYS_TOP" \
    -o "${DENSITY_DIR}/density_seg${seg}.tpr" \
    -po "${DENSITY_DIR}/mdout_density${seg}.mdp" \
    -maxwarn 5 \
    > "${DENSITY_DIR}/grompp_density${seg}.log" 2>&1

  $GMX mdrun \
    -deffnm "${DENSITY_DIR}/density_seg${seg}" \
    -x "${SCRATCH_DIR}/${OUTBASE}_density_seg${seg}.xtc" \
    -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu -pin on \
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
  LAST_CPT="${DENSITY_DIR}/density_seg${seg}.cpt"
done

if (( seg == DENSITY_MAX_SEG )); then
  echo "[WARNING] Reached DENSITY_MAX_SEG=${DENSITY_MAX_SEG} without convergence; using last segment"
fi

DENSITY_GRO="${DENSITY_DIR}/density_seg${seg}.gro"
DENSITY_CPT="${DENSITY_DIR}/density_seg${seg}.cpt"
[[ -f "$DENSITY_GRO" && -f "$DENSITY_CPT" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 6 — Relax (unrestrained NPT, production ensemble) — OPTIONAL
###############################################################################
# OFF by default (RELAX_NS=0). Restraints are released either HERE (RELAX_NS>0)
# or, when this step is skipped, at the start of production.
#
# Skipping it (default) makes the analyzed production trajectory CAPTURE the protein
# relaxing away from its restrained (designed) pose — the signal for stability /
# flexibility / variant-comparison. Enable it (RELAX_NS>0) when you instead want
# production to be a pre-equilibrated ensemble (e.g. bound-state sampling) and would
# rather drop the relaxation transient up front than trim it later. Either way
# production releases restraints (its mdp has no -DPOSRES) and keeps the same
# ensemble (V-rescale + C-rescale), so there is no thermostat/barostat switch.
CURRENT_STEP=6
if (( RELAX_NSTEPS > 0 )); then
echo "===== STEP 6: Relax (unrestrained NPT, ${RELAX_NS} ns) ====="

RELAX_DIR="${OUTDIR}/relax"
mkdir -p "$RELAX_DIR"
cat > "${RELAX_DIR}/relax.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${RELAX_NSTEPS}
nstxout-compressed = 5000
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
ref-t       = ${T_SIM} ${T_SIM}

pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 1.0
ref-p       = 1.0
compressibility = 4.5e-5

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = no
continuation = yes
EOF

$GMX grompp \
  -f "${RELAX_DIR}/relax.mdp" \
  -c "$DENSITY_GRO" \
  -t "$DENSITY_CPT" \
  -p "$SYS_TOP" \
  -o "${RELAX_DIR}/relax.tpr" \
  -po "${RELAX_DIR}/mdout_relax.mdp" \
  -maxwarn 5 \
  > "${RELAX_DIR}/grompp.log" 2>&1

$GMX mdrun \
  -deffnm "${RELAX_DIR}/relax" \
  -x "${SCRATCH_DIR}/relax.xtc" \
  -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu -pin on \
  > "${OUTDIR}/logs/mdrun_relax.log" 2>&1

RELAX_GRO="${RELAX_DIR}/relax.gro"
RELAX_CPT="${RELAX_DIR}/relax.cpt"
[[ -f "$RELAX_GRO" && -f "$RELAX_CPT" ]] || { echo "[ERROR] Relax (unrestrained NPT) failed"; exit 1; }
echo "[OK] Relax (unrestrained NPT) done"
else
echo "===== STEP 6: relax (unrestrained NPT) SKIPPED (RELAX_NS=${RELAX_NS}) ====="
echo "[INFO] Restraints release at production start; production will capture the"
echo "[INFO] relaxation from the restrained (designed) pose."
# Production continues straight from the last restrained NPT-density segment.
RELAX_GRO="$DENSITY_GRO"
RELAX_CPT="$DENSITY_CPT"
fi

###############################################################################
# STEP 7 — Production NPT (unrestrained, V-rescale + C-rescale)
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Production NPT ====="
echo "[INFO] ${TOTAL_NS} ns at ${T_SIM} K, frame every ${TRAJ_PS} ps"

PROD_DIR="${OUTDIR}/prod"
cat > "${PROD_DIR}/md.mdp" <<EOF
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
ref-t       = ${T_SIM} ${T_SIM}

pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 1.0
ref-p       = 1.0
compressibility = 4.5e-5

constraints = h-bonds
constraint-algorithm = LINCS

gen-vel     = no
continuation = yes

nstcomm     = 500
EOF

# Continue from the unrestrained eq if it ran, else directly from the last
# restrained NPT-density segment — in which case restraints release here, since
# md.mdp has no -DPOSRES (RELAX_GRO/RELAX_CPT are set to either source in STEP 6).
$GMX grompp \
  -f "${PROD_DIR}/md.mdp" \
  -c "$RELAX_GRO" \
  -t "$RELAX_CPT" \
  -p "$SYS_TOP" \
  -o "${PROD_DIR}/md.tpr" \
  -po "${PROD_DIR}/mdout_md.mdp" \
  -maxwarn 5 \
  > "${PROD_DIR}/grompp.log" 2>&1

$GMX mdrun \
  -deffnm "${PROD_DIR}/md" \
  -x "${SCRATCH_DIR}/md.xtc" \
  -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu -pin on \
  > "${OUTDIR}/logs/mdrun_md.log" 2>&1

[[ -f "${PROD_DIR}/md.gro" ]] || { echo "[ERROR] Production failed"; exit 1; }
echo "[OK] Production done"

###############################################################################
# STEP 8 — Finalize outputs & symlinks
###############################################################################
CURRENT_STEP=8
echo "===== STEP 8: Finalizing outputs ====="

# Symlink large trajectory files from scratch into the output tree
for f in "${SCRATCH_DIR}/heat.xtc" "${SCRATCH_DIR}/relax.xtc" "${SCRATCH_DIR}/md.xtc"; do
  [[ -f "$f" ]] && ln -sf "$f" "${OUTDIR}/trajectories/$(basename "$f")"
done
for f in "${SCRATCH_DIR}/${OUTBASE}_density_seg"*.xtc; do
  [[ -f "$f" ]] && ln -sf "$f" "${OUTDIR}/trajectories/$(basename "$f")"
done

# Export final PDB
echo 0 | $GMX trjconv \
  -s "${PROD_DIR}/md.tpr" \
  -f "${PROD_DIR}/md.gro" \
  -o "${OUTDIR}/${OUTBASE}_final.pdb" \
  -pbc whole \
  > /dev/null 2>&1 || true

echo "[OK] Outputs finalized"

###############################################################################
# STEP 9 — Write parameters log
###############################################################################
CURRENT_STEP=9
ELAPSED_S=$(( $(date +%s) - CURRENT_TIME ))

for dest in "${OUTDIR}/parameters.txt" "${SCRATCH_DIR}/parameters.txt"; do
  cat > "$dest" <<EOF
Elapsed seconds: $ELAPSED_S
Job Name: ${SLURM_JOB_NAME:-md}
Job ID: ${SLURM_JOB_ID:-N/A}
Total ns: $TOTAL_NS
Temperature: ${T_SIM} K
Force field: $FF
Water model: $WATER
Box shape: $BOX_SHAPE
Box buffer (nm): $BOX_BUFFER
Neutralize: $NEUTRALIZE
Salt molarity: $SALT_MOLAR
Cutoff (nm): $CUTOFF_NM
tau_t (ps): $TAU_T
Time step (ps): $DT_PS
Trajectory frame interval (ps): $TRAJ_PS
Production steps: $PROD_NSTEPS
Heat ns: $HEAT_NS
Density convergence segments: $seg
Relax (unrestrained NPT) ns: $RELAX_NS  (0 = disabled; restraints released at production start)
Production ensemble: NPT (V-rescale + C-rescale)
OMP threads: $NTOMP
Output directory: $OUTDIR
Scratch directory: $SCRATCH_DIR
GPUs: ${SLURM_GPUS_ON_NODE:-unknown}
EOF
done

###############################################################################
# STEP 10 — Post-analysis (PBC fix + strip/align, then RMSD/Rg/RMSF/DSSP + plots)
###############################################################################
CURRENT_STEP=10
echo "===== STEP 10: Post-analysis ====="

ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"

# The whole post-analysis lives in run_analysis.sh (auto-detects MD vs REMD from
# the job layout, activates the Python env, and echoes each [CMD] before running)
# so it can be re-run by hand later without resubmitting the simulation:
#   bash scripts/analysis/run_analysis.sh "$OUTDIR"
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" \
  || echo "[WARN] run_analysis.sh failed — re-run the command above"

###############################################################################
# STEP 11 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:    $OUTDIR"
echo "[INFO] Trajectory links:    ${OUTDIR}/trajectories/"
echo "[INFO] Actual trajectories: $SCRATCH_DIR"
echo "[INFO] Production log:       ${PROD_DIR}/md.log"
echo "[INFO] Processed traj:       ${OUTDIR}/analysis/md_stripped_aligned.xtc"
echo "[INFO] Analysis plots:       ${OUTDIR}/analysis/md_{rmsd,rg,rmsf,dssp}.png"
echo "=========================================="
