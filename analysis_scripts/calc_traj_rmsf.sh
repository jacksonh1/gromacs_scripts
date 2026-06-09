#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# calc_traj_rmsf.sh — per-residue backbone RMSF
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: works on any protein-only structure + trajectory pair,
# so it serves both the plain-MD and T-REMD pipelines unchanged.
#
# Usage:
#   bash calc_traj_rmsf.sh STRUCT TRAJ OUT
#
#   STRUCT  protein-only structure (.gro/.pdb/.tpr); topology for the analysis
#   TRAJ    protein-only trajectory (.xtc), already PBC-fixed and aligned
#   OUT     output .xvg  (per-residue RMSF)
#
# Computes the root-mean-square fluctuation of Backbone atoms, averaged per
# residue (-res), over the whole trajectory.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash calc_traj_rmsf.sh STRUCT TRAJ OUT}"
TRAJ="${2:?Usage: bash calc_traj_rmsf.sh STRUCT TRAJ OUT}"
OUT="${3:?Usage: bash calc_traj_rmsf.sh STRUCT TRAJ OUT}"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  set +u; source "$GMXRC"; set -u
fi

if command -v gmx_mpi &>/dev/null; then GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then GMX="gmx"
else echo "[ERROR] No GROMACS binary found. Source your GMXRC or set GROMACS_SCRIPTS_DIR."; exit 1
fi

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TRAJ" ]]   || { echo "[ERROR] Trajectory not found: $TRAJ"; exit 1; }
mkdir -p "$(dirname "$OUT")"

echo "[INFO] Per-residue RMSF: Backbone → $OUT"

# gmx rmsf prompts once for the analysis group; -res averages per residue.
printf "Backbone\n" | $GMX rmsf \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -o "$OUT" \
  -res

echo "[OK] RMSF written to: $OUT"
