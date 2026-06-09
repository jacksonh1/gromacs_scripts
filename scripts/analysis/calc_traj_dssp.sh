#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# calc_traj_dssp.sh — DSSP secondary structure over time
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: works on any protein-only structure + trajectory pair,
# so it serves both the plain-MD and T-REMD pipelines unchanged.
#
# Usage:
#   bash calc_traj_dssp.sh STRUCT TRAJ OUT
#
#   STRUCT  protein-only structure (.gro/.pdb/.tpr); topology for the analysis
#   TRAJ    protein-only trajectory (.xtc), already PBC-fixed and aligned
#   OUT     output .dat — one row per frame, one single-letter SS code per residue
#
# Uses gmx dssp (built into GROMACS 2023+, no external dssp binary needed).
# Hydrogens come from the structure (-hmode gromacs); the protein-only GRO from
# pdb2gmx already has them. Render with plot_dssp.py.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash calc_traj_dssp.sh STRUCT TRAJ OUT}"
TRAJ="${2:?Usage: bash calc_traj_dssp.sh STRUCT TRAJ OUT}"
OUT="${3:?Usage: bash calc_traj_dssp.sh STRUCT TRAJ OUT}"

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

echo "[INFO] DSSP secondary structure: Protein → $OUT"

# -sel takes a selection (not an interactive index group); Protein is the chain.
$GMX dssp \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -o "$OUT" \
  -sel "Protein"

echo "[OK] DSSP secondary structure written to: $OUT"
