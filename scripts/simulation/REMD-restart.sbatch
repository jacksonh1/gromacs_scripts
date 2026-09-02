#!/bin/bash
#SBATCH --job-name=remd-restart
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --gres=gpu:l40s:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# Resume a checkpointed T-REMD production run (companion to REMD-gromacs.sbatch)
# ---------------------------------------------------------------------------
# When a production REMD job stops early — SLURM TIMEOUT, a node failure, a
# manual scancel — every replica leaves a prod/repNNN/remd.cpt checkpoint. This
# script continues each replica from its checkpoint to the nsteps already baked
# into remd.tpr, appending to the existing scratch trajectory (one continuous
# run). The build/EM/density/equilibration stages are NOT repeated — the .tpr
# files and checkpoints from the original job carry everything they need.
#
# It does NOT extend the run beyond the original length. To run LONGER than the
# original TOTAL_NS, first `gmx convert-tpr -s remd.tpr -o remd.tpr -until <ps>`
# in each rep dir, then submit this script.
#
# Usage — match the ORIGINAL job's -n (replicas) and --gres:
#   sbatch -n 48 --gres=gpu:l40s:4 -t 4:00:00 \
#     --export=ALL,GROMACS_SCRIPTS_DIR=/abs/.../gromacs_REMD/scripts/simulation,\
# OUTDIR=/abs/path/to/<run>_<hash> \
#     /abs/.../scripts/simulation/REMD-restart.sbatch
#
# Or copy example/submit_jobs/submit_REMD_restart.sh, set OUTDIR, and run it.
#
# REPLEX_STEPS (the -replex flag) is not stored in the .tpr, so it is recovered
# from OUTDIR/logs/mdrun_remd.log. Override with --export ...,REPLEX_STEPS=500
# if that log is missing.
###############################################################################

set -euo pipefail
echo "===== REMD RESTART START $(date) ====="

# ── Required inputs ──────────────────────────────────────────────────────────
: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — the scripts/simulation dir of gromacs_REMD}"
: "${OUTDIR:?OUTDIR is not set — the existing run directory (contains prod/, logs/)}"
PROD_DIR="${OUTDIR}/prod"
[[ -d "$PROD_DIR" ]] || { echo "[ERROR] no prod/ under OUTDIR=$OUTDIR"; exit 1; }

# ── Environment (identical to REMD-gromacs.sbatch STEP 0) ────────────────────
source "${GROMACS_SCRIPTS_DIR}/../../site_config.sh"

module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"

set +u
source "$GMXRC"
set -u
GMX="${GMX:-gmx_mpi}"
command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }
echo "[OK] GROMACS found: $($GMX --version 2>&1 | head -1)"

# ── Discover replicas; every one must have a .tpr AND a checkpoint ───────────
mapfile -t REPDIRS < <(ls -d "${PROD_DIR}"/rep* 2>/dev/null | sort)
REPLICAS=${#REPDIRS[@]}
(( REPLICAS > 0 )) || { echo "[ERROR] no prod/rep* directories in $PROD_DIR"; exit 1; }

# One MPI rank per replica: SLURM_NTASKS must match, or mdrun -multidir mismatches.
if [[ -n "${SLURM_NTASKS:-}" && "$SLURM_NTASKS" -ne "$REPLICAS" ]]; then
  echo "[ERROR] SLURM_NTASKS=$SLURM_NTASKS but found $REPLICAS replicas — resubmit with -n $REPLICAS"
  exit 1
fi

for d in "${REPDIRS[@]}"; do
  [[ -f "$d/remd.tpr" ]] || { echo "[ERROR] missing $d/remd.tpr"; exit 1; }
  [[ -f "$d/remd.cpt" ]] || { echo "[ERROR] missing $d/remd.cpt — nothing to resume in $d"; exit 1; }
done

# ── Recover the exchange interval (-replex, in steps) ────────────────────────
# Not stored in the .tpr; read it from the original run's stdout log so the
# resumed exchange schedule matches. Override via the REPLEX_STEPS env var.
REPLEX_STEPS="${REPLEX_STEPS:-}"
if [[ -z "$REPLEX_STEPS" ]]; then
  ORIG_LOG="${OUTDIR}/logs/mdrun_remd.log"
  [[ -f "$ORIG_LOG" ]] || {
    echo "[ERROR] cannot auto-detect -replex: no $ORIG_LOG. Pass --export ...,REPLEX_STEPS=<steps>"; exit 1; }
  REPLEX_STEPS=$(grep -oE -- '-replex[[:space:]]+[0-9]+' "$ORIG_LOG" | head -1 | grep -oE '[0-9]+$' || true)
  [[ -n "$REPLEX_STEPS" ]] || {
    echo "[ERROR] no '-replex N' found in $ORIG_LOG. Pass --export ...,REPLEX_STEPS=<steps>"; exit 1; }
fi

# ── Recover the replica-exchange seed (-reseed) ──────────────────────────────
# The exchange RNG is counter-based on (seed, step), so a seeded original run only
# reproduces if the resume passes the SAME -reseed. Recover it from the original log;
# absent (older seedless run) or no log ⇒ -1 (random), matching that run's behaviour.
# Override via the SEED env var.
RESEED="${SEED:-}"
if [[ -z "$RESEED" && -f "${OUTDIR}/logs/mdrun_remd.log" ]]; then
  RESEED=$(grep -oE -- '-reseed[[:space:]]+-?[0-9]+' "${OUTDIR}/logs/mdrun_remd.log" | head -1 \
             | grep -oE -- '-?[0-9]+$' || true)
fi
RESEED="${RESEED:--1}"

echo "[OK] Resuming $REPLICAS replicas from checkpoint (-replex $REPLEX_STEPS, -reseed $RESEED)"

# ── Resume production ─────────────────────────────────────────────────────────
# -cpi remd.cpt enables continuation; with -multidir + -deffnm it resolves to
# each replica's own remd.cpt/.tpr/.xtc/.edr/.log. Append (the default) restores
# from the checkpoint and truncates the trajectory/energy/log back to that step
# before continuing. Under the folder-symlink model prod/ is a symlink into scratch
# (or a real dir when SYMLINK_BULK=0), so remd.xtc already lives at prod/rep*/remd.xtc.
# GPU/pin flags mirror REMD-gromacs.sbatch STEP 9.
RESTART_LOG="${OUTDIR}/logs/mdrun_remd_restart_${SLURM_JOB_ID:-manual}.log"

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir "${REPDIRS[@]}" \
  -deffnm remd \
  -cpi remd.cpt \
  -replex "$REPLEX_STEPS" \
  -reseed "$RESEED" \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "$RESTART_LOG" 2>&1

# ── Verify every replica reached the end (wrote its final structure) ──────────
for d in "${REPDIRS[@]}"; do
  [[ -f "$d/remd.gro" ]] || {
    echo "[ERROR] $d/remd.gro missing — production still incomplete (resubmit to continue)"; exit 1; }
done
echo "[OK] REMD production complete for all $REPLICAS replicas"

# ── Provenance: APPEND a restart record to the original parameters.txt ───────
# A restart can run under a different GROMACS build or a changed repo than the run
# that created this directory, so the record is appended rather than overwriting —
# parameters.txt ends up carrying the original run plus one entry per resume.
# finalize_outputs.sh copies parameters.txt into the scratch archive, so appending
# to OUTDIR is enough. Fail-soft probe (GOTCHAS: `VAR=$(cmd | grep ...)` + `set -e`).
set +e
GMX_VERSION="$($GMX --version 2>/dev/null | sed -n 's/^GROMACS version:[[:space:]]*//p' | head -1)"
SITE_CONFIG_PATH="$(cd "${GROMACS_SCRIPTS_DIR}/../.." && pwd)/site_config.sh"
set -e
GMX_VERSION="${GMX_VERSION:-unknown}"

cat >> "${OUTDIR}/parameters.txt" <<EOF

# ── Restart (job ${SLURM_JOB_ID:-manual}, $(date '+%Y-%m-%d %H:%M:%S')) ──
Engine: ${GROMACS_SCRIPTS_DIR}/REMD-restart.sbatch
Site config: $SITE_CONFIG_PATH
GROMACS: $GMX_VERSION
GMXRC: $GMXRC
GMXLIB: ${GMXLIB:-(unset)}
Replicas resumed: $REPLICAS
Exchange interval (steps): $REPLEX_STEPS
RNG seed (-reseed): $RESEED
Restart mdrun log: $RESTART_LOG
Submitted from: ${SLURM_SUBMIT_DIR:-N/A}
EOF
echo "[OK] Restart provenance appended to ${OUTDIR}/parameters.txt"

# ── Finalize: export final PDB (mirrors STEP 10) ─────────────────────────────
# No trajectory re-linking: under the folder-symlink model remd.xtc is already at
# prod/rep*/remd.xtc (real on scratch via the prod/ symlink, or real in OUTDIR).
echo 0 | $GMX trjconv \
  -s "${PROD_DIR}/rep000/remd.tpr" \
  -f "${PROD_DIR}/rep000/remd.gro" \
  -o "${OUTDIR}/$(basename "$OUTDIR")_final_rep000.pdb" \
  -pbc whole \
  > /dev/null 2>&1 || true

# ── Post-analysis (reuse the toolset entry point; do not duplicate its logic) ─
ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR 000"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" 000 \
  || echo "[WARN] run_analysis.sh failed — re-run the command above"

# ── Refresh the scratch archive (if the run was offloaded to scratch) ─────────
# prod/ is a symlink into scratch when SYMLINK_BULK=1; recover SCRATCH_DIR as its
# parent and re-copy the small OUTDIR dirs so the archive reflects the resumed run.
# When prod/ is a real dir (SYMLINK_BULK=0), there is no scratch archive → skip.
if [[ -L "${OUTDIR}/prod" ]]; then
  SCRATCH_DIR="$(dirname "$(readlink -f "${OUTDIR}/prod")")"
  echo "[CMD] bash ${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh $OUTDIR $SCRATCH_DIR"
  bash "${GROMACS_SCRIPTS_DIR}/finalize_outputs.sh" "$OUTDIR" "$SCRATCH_DIR" \
    || echo "[WARN] finalize_outputs.sh failed — scratch archive may be incomplete; re-run the command above"
else
  SCRATCH_DIR=""
fi

echo ""
echo "=========================================="
echo "===== RESTART COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:  $OUTDIR"
echo "[INFO] Restart mdrun log: $RESTART_LOG"
echo "[INFO] Exchange log:      ${PROD_DIR}/rep000/remd.log"
[[ -n "$SCRATCH_DIR" ]] && echo "[INFO] Scratch archive:   $SCRATCH_DIR"
echo "[INFO] Processed traj:    ${OUTDIR}/analysis/remd_rep000_stripped_aligned.xtc"
echo "=========================================="
