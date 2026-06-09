#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_chain_rmsd.sh — backbone RMSD of ONE chain, fit to itself
# ─────────────────────────────────────────────────────────────────────────────
# Per-chain counterpart of calc_traj_rmsd.sh, for the multi-chain pipeline. Fits
# and computes RMSD on a single chain's backbone (group from the chain index), so
# the result is that chain's internal drift — independent of how the chains sit
# relative to each other. This is immune to the inter-chain PBC ambiguity and to
# (un)binding: it answers "does this chain keep its fold?".
#
# Usage:
#   bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT
#
#   STRUCT  protein-only reference (.gro); RMSD reference and topology
#   TRAJ    protein-only trajectory (.xtc), PBC-fixed and aligned
#   NDX     index file with the chain backbone group (from multichain_chain_index.py)
#   GROUP   chain backbone group name, e.g. ChainA_Backbone
#   OUT     output .xvg
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT}"
TRAJ="${2:?Usage: bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT}"
NDX="${3:?Usage: bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT}"
GROUP="${4:?Usage: bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT}"
OUT="${5:?Usage: bash multichain_chain_rmsd.sh STRUCT TRAJ NDX GROUP OUT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then source "$SITE_CONFIG"; set +u; source "$GMXRC"; set -u; fi
if command -v gmx_mpi &>/dev/null; then GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then GMX="gmx"
else echo "[ERROR] No GROMACS binary found."; exit 1; fi

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TRAJ" ]]   || { echo "[ERROR] Trajectory not found: $TRAJ"; exit 1; }
[[ -f "$NDX" ]]    || { echo "[ERROR] Index not found: $NDX"; exit 1; }
mkdir -p "$(dirname "$OUT")"

echo "[INFO] Chain RMSD: fit ${GROUP}, compute ${GROUP} → $OUT"
# gmx rms prompts twice: (1) least-squares fit group, (2) RMSD group.
printf "%s\n%s\n" "$GROUP" "$GROUP" | $GMX rms \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -n "$NDX" \
  -o "$OUT"

echo "[OK] Chain RMSD written to: $OUT"
