#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# strip_and_align_trajectory.sh — strip waters/ions and align protein
# ─────────────────────────────────────────────────────────────────────────────
# Requires process_trajectory.sh to have been run first (needs the
# PBC-corrected trajectory in OUTDIR/analysis/).
#
# Usage:
#   bash strip_and_align_trajectory.sh OUTDIR [REP] [REF_GRO]
#
#   OUTDIR    path to the job output directory
#   REP       replica index (default: 000 — the 300 K replica)
#   REF_GRO   reference structure for alignment (default: first frame of the
#             PBC-corrected trajectory)
#
# Output written to: OUTDIR/analysis/
#   remd_rep<REP>_ref_full.gro — full-system first frame (internal; used as -s
#                                for the fitting call to match XTC atom count)
#   remd_rep<REP>_ref.gro      — protein-only first frame (use this as the
#                                reference structure for RMSD, clustering, etc.)
#   remd_rep<REP>_fit.xtc      — protein-only trajectory, backbone-aligned
#
# Alignment: backbone fit (N, CA, C, O) to REF_GRO; protein atoms output.
# Stripping: waters and ions are dropped by selecting the Protein output group.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

OUTDIR="${1:?Usage: bash strip_and_align_trajectory.sh OUTDIR [REP] [REF_GRO]}"
REP="${2:-000}"
REF_GRO="${3:-}"

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

# ── Input files ───────────────────────────────────────────────────────────────
TPR="${OUTDIR}/prod/rep${REP}/remd.tpr"
PBC_XTC="${OUTDIR}/analysis/remd_rep${REP}_pbc.xtc"

[[ -f "$TPR" ]]     || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -f "$PBC_XTC" ]] || { echo "[ERROR] PBC-corrected XTC not found: $PBC_XTC"
                          echo "        Run process_trajectory.sh first."; exit 1; }

ANALYSIS_DIR="${OUTDIR}/analysis"
# Full-system GRO — used as -s for gmx trjconv fitting (must match XTC atom count)
REF_FULL="${ANALYSIS_DIR}/remd_rep${REP}_ref_full.gro"
# Protein-only GRO — the reference structure for downstream analysis (RMSD, etc.)
REF_PROTEIN="${ANALYSIS_DIR}/remd_rep${REP}_ref.gro"
OUT_XTC="${ANALYSIS_DIR}/remd_rep${REP}_fit.xtc"

# ── Reference structure ───────────────────────────────────────────────────────
# The -s flag passed to gmx trjconv must have the SAME atom count as the XTC
# (full system: protein + water + ions). Passing a protein-only GRO silently
# truncates the trajectory read and causes the Jacobi rotation fitting to fail.
#
# We therefore extract two versions of frame 0:
#   _ref_full.gro  — full system, used internally as -s for the fitting call
#   _ref.gro       — protein-only, the reference for downstream analysis
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

echo ""
echo "[OK] Done."
echo "  Reference (protein-only): $REF_PROTEIN"
echo "  Trajectory (protein-only, aligned): $OUT_XTC"
