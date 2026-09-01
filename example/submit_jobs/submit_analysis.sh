#!/usr/bin/env bash
# Re-run the post-analysis on an EXISTING job directory (no simulation is re-run).
# Copy this file, set OUTDIR, and run it.
# Usage: bash submit_analysis.sh
#
# Use this when an analysis script changed or was added, when the job's analysis
# step failed but the trajectory is fine, or to re-run with different knobs.
# MD vs T-REMD vs REST2 is auto-detected from the job layout — one script for all
# three. The trajectory must still exist: under SYMLINK_BULK=1 prod/ is a symlink
# into scratch, so a purged scratch dir means there is nothing left to analyse.
#
# On a login node or an interactive allocation you can skip SLURM entirely:
#   bash ../../scripts/analysis/run_analysis.sh "$OUTDIR"

export GROMACS_SCRIPTS_DIR="/orcd/pool/004/jhalpin/09-fragfold/RELE_simulations/gromacs_REMD/scripts/simulation"

# === Job to analyse ===
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_MD/helix_fusion-2ns-MD-300K-NPT"

# === Replica slot (T-REMD / REST2 only; ignored for plain MD) ===
# 000 is the lowest temperature (T-REMD) / lambda=1 (REST2) — the physical ensemble.
REP=000

# === Analysis knobs (optional) ===
# CLUSTER_CUTOFF  backbone-RMSD cutoff for conformational clustering, nm (default 0.20 = 2.0 Å)
# SHELL_NM        solvated-snapshot solvent shell, nm (default 0.5 = 5 Å) — MD only
# N_SNAPSHOTS     how many evenly-spaced solvated snapshots (default 5)   — MD only
CLUSTER_CUTOFF=0.20
SHELL_NM=0.5
N_SNAPSHOTS=5


sbatch \
  --export=ALL,GROMACS_SCRIPTS_DIR="$GROMACS_SCRIPTS_DIR",OUTDIR="$OUTDIR",REP="$REP",CLUSTER_CUTOFF="$CLUSTER_CUTOFF",SHELL_NM="$SHELL_NM",N_SNAPSHOTS="$N_SNAPSHOTS" \
  "${GROMACS_SCRIPTS_DIR}/../analysis/analysis.sbatch"
