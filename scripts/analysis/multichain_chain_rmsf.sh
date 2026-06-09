#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_chain_rmsf.sh — per-residue backbone RMSF of ONE chain
# ─────────────────────────────────────────────────────────────────────────────
# Per-chain counterpart of calc_traj_rmsf.sh, for the multi-chain pipeline.
# Computes per-residue RMSF over a single chain's backbone (group from the chain
# index), so per-chain flexibility is reported separately rather than concatenated.
#
# Usage:
#   bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT
#
#   STRUCT  protein-only structure (.gro); topology / averaging reference
#   TRAJ    protein-only trajectory (.xtc), PBC-fixed and aligned
#   NDX     index file with the chain backbone group (from multichain_chain_index.py)
#   GROUP   chain backbone group name, e.g. ChainA_Backbone
#   OUT     output .xvg
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT}"
TRAJ="${2:?Usage: bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT}"
NDX="${3:?Usage: bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT}"
GROUP="${4:?Usage: bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT}"
OUT="${5:?Usage: bash multichain_chain_rmsf.sh STRUCT TRAJ NDX GROUP OUT}"

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

echo "[INFO] Chain RMSF: per-residue over ${GROUP} → $OUT"
printf "%s\n" "$GROUP" | $GMX rmsf \
  -s "$STRUCT" \
  -f "$TRAJ" \
  -n "$NDX" \
  -res \
  -o "$OUT"

echo "[OK] Chain RMSF written to: $OUT"
