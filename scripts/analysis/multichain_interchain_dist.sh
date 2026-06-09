#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_interchain_dist.sh — minimum-image distance between two chains
# ─────────────────────────────────────────────────────────────────────────────
# The binding observable for a complex: the minimum distance between two chains
# over the trajectory (`gmx mindist`, which is minimum-image-aware, so it is the
# real contact distance regardless of PBC). A near-constant small value (~0.15 nm
# = contact) means the complex stayed bound; a rise toward ~half the box flags
# (un)binding AND that the box is too small for that separation.
#
# Usage:
#   bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT
#
#   STRUCT        protein-only structure (.gro)
#   TRAJ          protein-only trajectory (.xtc)
#   NDX           index with the two chain groups (from multichain_chain_index.py)
#   GROUP1 GROUP2 chain group names, e.g. ChainA ChainB
#   OUT           output .xvg (minimum distance vs time)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"
TRAJ="${2:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"
NDX="${3:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"
GROUP1="${4:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"
GROUP2="${5:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"
OUT="${6:?Usage: bash multichain_interchain_dist.sh STRUCT TRAJ NDX GROUP1 GROUP2 OUT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then source "$SITE_CONFIG"; set +u; source "$GMXRC"; set -u; fi
if command -v gmx_mpi &>/dev/null; then GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then GMX="gmx"
else echo "[ERROR] No GROMACS binary (gmx_mpi/gmx) on PATH. Source your GROMACS GMXRC or load the GROMACS module first."; exit 1; fi

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TRAJ" ]]   || { echo "[ERROR] Trajectory not found: $TRAJ"; exit 1; }
[[ -f "$NDX" ]]    || { echo "[ERROR] Index not found: $NDX"; exit 1; }
mkdir -p "$(dirname "$OUT")"

echo "[INFO] Inter-chain min distance: ${GROUP1} <-> ${GROUP2} → $OUT"
# gmx mindist prompts for the two groups; -od writes minimum distance vs time.
printf "%s\n%s\n" "$GROUP1" "$GROUP2" | $GMX mindist \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -n "$NDX" \
  -od "$OUT"

echo "[OK] Inter-chain distance written to: $OUT"
