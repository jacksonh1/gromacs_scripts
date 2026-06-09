#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# calc_traj_rg.sh — radius of gyration vs time
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: works on any protein-only structure + trajectory pair,
# so it serves both the plain-MD and T-REMD pipelines unchanged.
#
# Usage:
#   bash calc_traj_rg.sh STRUCT TRAJ OUT
#
#   STRUCT  protein-only structure (.gro/.pdb/.tpr); topology for the analysis
#   TRAJ    protein-only trajectory (.xtc), already PBC-fixed and aligned
#   OUT     output .xvg
#
# Computes the radius of gyration over the Protein group for each frame.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash calc_traj_rg.sh STRUCT TRAJ OUT}"
TRAJ="${2:?Usage: bash calc_traj_rg.sh STRUCT TRAJ OUT}"
OUT="${3:?Usage: bash calc_traj_rg.sh STRUCT TRAJ OUT}"

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

echo "[INFO] Radius of gyration: Protein → $OUT"

# gmx gyrate prompts once for the analysis group.
printf "Protein\n" | $GMX gyrate \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -o "$OUT"

echo "[OK] Radius of gyration written to: $OUT"
