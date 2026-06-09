#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# fix_PBC_strip_align.sh — PBC fix + strip waters/ions + alignment, no temp file
# ─────────────────────────────────────────────────────────────────────────────
# Orchestrator: calls fix_PBC.sh then strip_and_align_trajectory.sh, then removes
# the intermediate full-system PBC trajectory (<prefix>_pbc.xtc) to save disk.
#
# This is the ONE place that encodes the <prefix>_pbc.xtc intermediate naming
# convention; the base scripts it calls are layout-blind (explicit paths only).
#
# Use this in automated pipelines. Use fix_PBC.sh directly when you need to keep
# the full-system (solvent-included) PBC-corrected trajectory.
#
# Usage:
#   bash fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]
#
#   TPR         run input matching the XTC atom count (full system)
#   XTC         raw trajectory to process
#   OUT_PREFIX  output path stem (e.g. OUTDIR/analysis/md  or
#                                       OUTDIR/analysis/remd_rep000)
#   REF_GRO     optional full-system reference structure for alignment
#               (default: first frame of the trajectory)
#
# Output written to:
#   <OUT_PREFIX>_stripped_aligned.xtc   — protein-only trajectory, backbone-aligned (kept)
#   <OUT_PREFIX>_stripped_aligned.gro   — protein-only reference structure          (kept)
#   <OUT_PREFIX>_pbc.xtc   — full-system PBC intermediate              (deleted)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

TPR="${1:?Usage: bash fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
XTC="${2:?Usage: bash fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
OUT_PREFIX="${3:?Usage: bash fix_PBC_strip_align.sh TPR XTC OUT_PREFIX [REF_GRO]}"
REF_GRO="${4:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PBC_XTC="${OUT_PREFIX}_pbc.xtc"

bash "${SCRIPT_DIR}/fix_PBC.sh" "$TPR" "$XTC" "$PBC_XTC"
bash "${SCRIPT_DIR}/strip_and_align_trajectory.sh" "$TPR" "$PBC_XTC" "$OUT_PREFIX" ${REF_GRO:+"$REF_GRO"}

echo "[INFO] Removing intermediate PBC trajectory: $PBC_XTC"
rm -f "$PBC_XTC"

echo "[OK] Done. Intermediate _pbc.xtc removed."
echo "  Reference (protein-only): ${OUT_PREFIX}_stripped_aligned.gro"
echo "  Trajectory (protein-only, aligned): ${OUT_PREFIX}_stripped_aligned.xtc"
