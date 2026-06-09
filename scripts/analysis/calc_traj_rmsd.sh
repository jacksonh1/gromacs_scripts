#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# calc_traj_rmsd.sh — backbone RMSD vs time
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: works on any protein-only structure + trajectory pair,
# so it serves both the plain-MD and T-REMD pipelines unchanged.
#
# Usage:
#   bash calc_traj_rmsd.sh STRUCT TRAJ OUT
#
#   STRUCT  protein-only reference structure (.gro/.pdb/.tpr); used as both the
#           topology and the RMSD reference. Must match TRAJ's atom count.
#   TRAJ    protein-only trajectory (.xtc), already PBC-fixed and aligned
#   OUT     output .xvg
#
# Computes least-squares fit on Backbone and RMSD on Backbone (N, CA, C, O).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash calc_traj_rmsd.sh STRUCT TRAJ OUT}"
TRAJ="${2:?Usage: bash calc_traj_rmsd.sh STRUCT TRAJ OUT}"
OUT="${3:?Usage: bash calc_traj_rmsd.sh STRUCT TRAJ OUT}"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  set +u; source "$GMXRC"; set -u
fi

if command -v gmx_mpi &>/dev/null; then GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then GMX="gmx"
else echo "[ERROR] No GROMACS binary (gmx_mpi/gmx) on PATH. Source your GROMACS GMXRC or load the GROMACS module first."; exit 1
fi

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TRAJ" ]]   || { echo "[ERROR] Trajectory not found: $TRAJ"; exit 1; }
mkdir -p "$(dirname "$OUT")"

echo "[INFO] RMSD: fit Backbone, compute Backbone → $OUT"

# gmx rms prompts twice: (1) least-squares fit group, (2) RMSD group.
printf "Backbone\nBackbone\n" | $GMX rms \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -o "$OUT"

echo "[OK] RMSD written to: $OUT"
