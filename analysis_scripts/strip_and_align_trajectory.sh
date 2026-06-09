#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# strip_and_align_trajectory.sh — strip waters/ions and align protein
# ─────────────────────────────────────────────────────────────────────────────
# Layout-blind: takes explicit input paths and an output prefix, so the same
# script works for any pipeline (plain MD, T-REMD, …).
#
# Requires a PBC-corrected, full-system trajectory (run fix_PBC.sh first).
#
# Usage:
#   bash strip_and_align_trajectory.sh TPR PBC_XTC OUT_PREFIX [REF_GRO]
#
#   TPR         run input matching the PBC_XTC atom count (full system)
#   PBC_XTC     PBC-corrected, full-system trajectory (from fix_PBC.sh)
#   OUT_PREFIX  output path stem; the script writes a matched pair (same stem):
#                 <OUT_PREFIX>_stripped_aligned.xtc  protein-only, backbone-aligned
#                                                    trajectory (kept)
#                 <OUT_PREFIX>_stripped_aligned.gro  its first frame — the protein-only
#                                                    structure/reference for that
#                                                    trajectory's analysis (kept)
#               and one internal throwaway:
#                 <OUT_PREFIX>_frame0_fullsys.gro    full-system first frame, used as
#                                                    -s to match the XTC atom count,
#                                                    then deleted
#   REF_GRO     optional full-system reference structure to align to instead of
#               the first trajectory frame
#
# Alignment: backbone fit (N, CA, C, O) to the reference; protein atoms output.
# Stripping: waters and ions are dropped by selecting the Protein output group.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

TPR="${1:?Usage: bash strip_and_align_trajectory.sh TPR PBC_XTC OUT_PREFIX [REF_GRO]}"
PBC_XTC="${2:?Usage: bash strip_and_align_trajectory.sh TPR PBC_XTC OUT_PREFIX [REF_GRO]}"
OUT_PREFIX="${3:?Usage: bash strip_and_align_trajectory.sh TPR PBC_XTC OUT_PREFIX [REF_GRO]}"
REF_GRO="${4:-}"

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
[[ -f "$TPR" ]]     || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -f "$PBC_XTC" ]] || { echo "[ERROR] PBC-corrected XTC not found: $PBC_XTC"
                          echo "        Run fix_PBC.sh first."; exit 1; }

mkdir -p "$(dirname "$OUT_PREFIX")"

# Full-system GRO — used as -s for gmx trjconv fitting (must match XTC atom count)
REF_FULL="${OUT_PREFIX}_frame0_fullsys.gro"
# Protein-only GRO — the reference structure for downstream analysis (RMSD, etc.)
REF_PROTEIN="${OUT_PREFIX}_stripped_aligned.gro"
OUT_XTC="${OUT_PREFIX}_stripped_aligned.xtc"

# ── Reference structure ───────────────────────────────────────────────────────
# The -s flag passed to gmx trjconv must have the SAME atom count as the XTC
# (full system: protein + water + ions). Passing a protein-only GRO silently
# truncates the trajectory read and causes the Jacobi rotation fitting to fail.
#
# We therefore extract two versions of frame 0:
#   _frame0_fullsys.gro    — full system, used internally as -s for the fitting call
#   _stripped_aligned.gro  — protein-only, the reference for downstream analysis
#
# Custom: provide a full-system GRO/PDB as REF_GRO to align to a specific
# structure instead of the first trajectory frame.
if [[ -z "$REF_GRO" ]]; then
  echo "[INFO] Extracting first frame as alignment reference..."
  printf "System\n" | $GMX trjconv \
    -s "$TPR" -f "$PBC_XTC" -o "$REF_FULL" -dump 0
  printf "Protein\n" | $GMX trjconv \
    -s "$TPR" -f "$PBC_XTC" -o "$REF_PROTEIN" -dump 0
  REF_GRO="$REF_FULL"
  echo "[INFO] Reference (full system): $REF_FULL"
  echo "[INFO] Reference (protein):     $REF_PROTEIN"
else
  echo "[INFO] Using provided reference: $REF_GRO"
  [[ -f "$REF_GRO" ]] || { echo "[ERROR] Reference not found: $REF_GRO"; exit 1; }
  # Still need a protein-only reference for downstream analysis.
  printf "Protein\n" | $GMX trjconv \
    -s "$TPR" -f "$PBC_XTC" -o "$REF_PROTEIN" -dump 0
fi

# ── Strip and align ───────────────────────────────────────────────────────────
echo "[INFO] Fitting group: Backbone | Output group: Protein"
echo "[INFO] Output XTC: $OUT_XTC"
echo ""

# Two selections prompted in order:
#   1. Fitting group  → Backbone (N, CA, C, O)
#   2. Output group   → Protein
printf "Backbone\nProtein\n" | $GMX trjconv \
  -s "$REF_GRO" \
  -f "$PBC_XTC" \
  -o "$OUT_XTC" \
  -fit rot+trans

# The full-system reference is only needed as the -s for the fit call above.
rm -f "$REF_FULL"

echo ""
echo "[OK] Done."
echo "  Reference (protein-only): $REF_PROTEIN"
echo "  Trajectory (protein-only, aligned): $OUT_XTC"
