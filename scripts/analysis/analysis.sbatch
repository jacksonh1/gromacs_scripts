#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# Re-run the post-analysis on an EXISTING job directory
# ---------------------------------------------------------------------------
# The simulation engines already run the analysis as their last step, so this is
# for the cases where you want it again without resubmitting the simulation:
#   - an analysis script was fixed or added (e.g. new solvated snapshots)
#   - the job's analysis step failed but the trajectory is fine
#   - you want to re-run with different knobs (CLUSTER_CUTOFF, SHELL_NM, ...)
#   - you want a different T-REMD/REST2 replica slot (REP)
#
# No GPU is requested: every analysis step is CPU-only (gmx trjconv/rms/gyrate/
# rmsf/dssp, sklearn, matplotlib). On pi_keating that lands on the CPU-only nodes,
# so this never sits in the GPU queue behind a production run.
#
# All the actual work is run_analysis.sh — this file only supplies the SLURM
# allocation and the module/GROMACS environment a compute node needs. On a login
# node or an interactive allocation, skip sbatch entirely:
#
#   bash scripts/analysis/run_analysis.sh OUTDIR [REP]
#
# Intended usage:
#   Copy example/submit_jobs/submit_analysis.sh, set OUTDIR, run it.
#
# Advanced / direct usage:
#   sbatch --export=ALL,GROMACS_SCRIPTS_DIR=/path/scripts/simulation,OUTDIR=/path/job \
#     scripts/analysis/analysis.sbatch
###############################################################################

set -euo pipefail

echo "===== ANALYSIS JOB START $(date) ====="

# NOTE: SLURM copies an sbatch script to a temp path, so BASH_SOURCE[0] does NOT
# point into the repo. The repo location arrives as GROMACS_SCRIPTS_DIR — the same
# variable every submit script in this repo already sets.
: "${GROMACS_SCRIPTS_DIR:?GROMACS_SCRIPTS_DIR is not set — submit via example/submit_jobs/submit_analysis.sh}"
ANALYSIS_SCRIPTS_DIR="${GROMACS_SCRIPTS_DIR}/../analysis"

OUTDIR="${OUTDIR:?OUTDIR is not set — point it at a finished job directory}"
REP="${REP:-000}"          # T-REMD / REST2 replica slot; ignored for plain MD

[[ -d "$OUTDIR" ]] || { echo "[ERROR] Not a directory: $OUTDIR"; exit 1; }

source "${GROMACS_SCRIPTS_DIR}/../../site_config.sh"

module purge
module load gcc
module load "$CUDA_MODULE"
module load "$OPENMPI_MODULE"

# gmx_mpi links libhcoll.so.1 / libocoms.so.0, which some compute nodes lack under
# /opt/mellanox — without this it dies instantly with "cannot open shared object".
# Harmless no-op on nodes that already have them (see CLAUDE.md gotchas).
export LD_LIBRARY_PATH="${HCOLL_COMPAT_DIR}:${LD_LIBRARY_PATH:-}"

echo "[INFO] OUTDIR : $OUTDIR"
echo "[INFO] REP    : $REP  (T-REMD/REST2 only)"
echo "[INFO] Knobs  : CLUSTER_CUTOFF=${CLUSTER_CUTOFF:-0.20} SHELL_NM=${SHELL_NM:-0.5} N_SNAPSHOTS=${N_SNAPSHOTS:-5}"
echo ""

# run_analysis.sh sources GMXRC and activates the conda env itself, and auto-detects
# MD vs T-REMD vs REST2 from the job layout — nothing engine-specific belongs here.
echo "[CMD] bash ${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh $OUTDIR $REP"
bash "${ANALYSIS_SCRIPTS_DIR}/run_analysis.sh" "$OUTDIR" "$REP"

echo ""
echo "===== ANALYSIS JOB COMPLETED $(date) ====="
# NOT refreshed here: under SYMLINK_BULK=1 the original job copied analysis/ into the
# scratch archive via finalize_outputs.sh, and that copy is now stale. Refreshing it is
# a deliberate, separate step (the archive is a snapshot, and re-running it silently
# would overwrite whatever is there) — the scratch path is in OUTDIR/parameters.txt:
#   bash scripts/simulation/finalize_outputs.sh "$OUTDIR" "$SCRATCH_DIR"
echo "[INFO] Plots / metrics:    ${OUTDIR}/analysis/"
echo "[INFO] Clustering:         ${OUTDIR}/analysis/clustering/"
# `if`, not `[[ ... ]] &&` — under `set -e` a false test as the LAST command would
# exit nonzero and SLURM would report the whole job as FAILED.
if [[ -d "${OUTDIR}/solvated_snapshots" ]]; then
  echo "[INFO] Solvated snapshots: ${OUTDIR}/solvated_snapshots/"
fi
