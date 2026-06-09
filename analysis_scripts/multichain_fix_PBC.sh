#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_fix_PBC.sh — PBC fix for a MULTI-CHAIN complex (keeps chains together)
# ─────────────────────────────────────────────────────────────────────────────
# Multi-chain counterpart of fix_PBC.sh. For a complex, plain `-pbc mol` wraps
# each chain's centre of mass into the box *independently*, so the chains can land
# in different periodic images and the complex is split across a box boundary.
# This script keeps the chains in ONE image, in three per-frame passes:
#   1. -pbc whole    : make every molecule whole.
#   2. -pbc cluster  : pull the protein chains into the same periodic image
#                      (cluster group = Protein). NOTE: "-pbc cluster" is a
#                      PERIODIC-IMAGE operation — it is NOT conformational
#                      clustering (gmx cluster); this script handles PBC only.
#   3. -pbc mol -center -ur compact : centre the protein, compact box for viewing.
# All three are per-frame (no frame-to-frame comparison), so this is REMD-safe.
#
# ASSUMES the complex stays within ~half the (minimum) box vector. If chains
# dissociate further, the periodic image is genuinely ambiguous (and the box is
# too small — a minimum-image violation). run_analysis emits an inter-chain
# minimum-distance curve so that case is visible.
#
# Layout-blind: explicit input/output paths (same contract as fix_PBC.sh).
#
# Usage:
#   bash multichain_fix_PBC.sh TPR XTC OUT_PBC_XTC
#
#   TPR          run input matching the XTC atom count (full system)
#   XTC          raw trajectory to correct
#   OUT_PBC_XTC  output path for the PBC-corrected (full-system) trajectory
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

TPR="${1:?Usage: bash multichain_fix_PBC.sh TPR XTC OUT_PBC_XTC}"
XTC="${2:?Usage: bash multichain_fix_PBC.sh TPR XTC OUT_PBC_XTC}"
OUT_XTC="${3:?Usage: bash multichain_fix_PBC.sh TPR XTC OUT_PBC_XTC}"

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

# ── Validate inputs ───────────────────────────────────────────────────────────
[[ -f "$TPR" ]] || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -f "$XTC" ]] || { echo "[ERROR] XTC not found: $XTC"; exit 1; }
mkdir -p "$(dirname "$OUT_XTC")"

echo "[INFO] Input TPR: $TPR"
echo "[INFO] Input XTC: $XTC"
echo "[INFO] Output:    $OUT_XTC"
echo ""

# ── PBC correction (3 per-frame passes; see header) ───────────────────────────
# Intermediates derived from the output name; cleaned up on exit (even on error).
TMP_WHOLE="${OUT_XTC%.xtc}.tmp_whole.xtc"
TMP_CLUST="${OUT_XTC%.xtc}.tmp_clust.xtc"
trap 'rm -f "$TMP_WHOLE" "$TMP_CLUST"' EXIT

# 1. make molecules whole (output group: System)
printf "System\n" | $GMX trjconv -s "$TPR" -f "$XTC" -o "$TMP_WHOLE" -pbc whole

# 2. bring the protein chains into one image (cluster group: Protein, output: System)
printf "Protein\nSystem\n" | $GMX trjconv -s "$TPR" -f "$TMP_WHOLE" -o "$TMP_CLUST" -pbc cluster

# 3. centre protein, compact box (centre group: Protein, output: System)
printf "Protein\nSystem\n" | $GMX trjconv -s "$TPR" -f "$TMP_CLUST" -o "$OUT_XTC" \
  -pbc mol -center -ur compact

echo ""
echo "[OK] PBC-corrected (chains kept together) trajectory written to: $OUT_XTC"
