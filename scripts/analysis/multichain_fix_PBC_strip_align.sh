#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_fix_PBC_strip_align.sh — multi-chain PBC fix + strip + align
# ─────────────────────────────────────────────────────────────────────────────
# Multi-chain counterpart of fix_PBC_strip_align.sh. Orchestrator: calls
# multichain_fix_PBC.sh (keeps the complex's chains in one periodic image) then the
# shared, chain-agnostic strip_and_align_trajectory.sh, then removes the
# intermediate full-system PBC trajectory.
#
# The strip/align step is identical to the single-chain path — it only does a
# backbone least-squares fit and frame dumps (no PBC), so it preserves the
# already-clustered complex.
#
# Usage:
#   bash multichain_fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]
#
# Output (same names as the single-chain path, so analysis mirrors it):
#   <OUT_PREFIX>_stripped_aligned.xtc   protein-only, backbone-aligned (kept)
#   <OUT_PREFIX>_stripped_aligned.gro   protein-only reference          (kept)
#   <OUT_PREFIX>_pbc.xtc                full-system PBC intermediate    (deleted)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

TPR="${1:?Usage: bash multichain_fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
XTC="${2:?Usage: bash multichain_fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
OUT_PREFIX="${3:?Usage: bash multichain_fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
REF_GRO="${4:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PBC_XTC="${OUT_PREFIX}_pbc.xtc"

bash "${SCRIPT_DIR}/multichain_fix_PBC.sh" "$TPR" "$XTC" "$PBC_XTC"
bash "${SCRIPT_DIR}/strip_and_align_trajectory.sh" "$TPR" "$PBC_XTC" "$OUT_PREFIX" ${REF_GRO:+"$REF_GRO"}

echo "[INFO] Removing intermediate PBC trajectory: $PBC_XTC"
rm -f "$PBC_XTC"

echo "[OK] Done. Intermediate _pbc.xtc removed."
echo "  Reference (protein-only): ${OUT_PREFIX}_stripped_aligned.gro"
echo "  Trajectory (protein-only, aligned): ${OUT_PREFIX}_stripped_aligned.xtc"
