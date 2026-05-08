#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# process_trajectory.sh — fix PBC artifacts in the 300 K REMD trajectory
# ─────────────────────────────────────────────────────────────────────────────
# Usage:
#   bash process_trajectory.sh OUTDIR [REP]
#
#   OUTDIR   path to the job output directory (contains prod/, trajectories/)
#   REP      replica index to process (default: 000 — the 300 K replica)
#
# Output written to: OUTDIR/analysis/remd_rep<REP>_pbc.xtc
#
# PBC treatment: -pbc mol -center -ur compact
#   - mol:     makes each molecule whole and puts its COM inside the box.
#              Per-frame operation — safe for REMD trajectories.
#   - center:  centers the protein in the box.
#   - compact: renders the dodecahedron box as a compact shape; without this
#              the box looks "exploded" in VMD/PyMOL.
#
# WARNING: -pbc nojump is intentionally NOT used here. nojump compares
# consecutive frames to detect box-boundary crossings, but in T-REMD,
# coordinate exchanges cause discontinuous jumps between frames that nojump
# would misinterpret and incorrectly "fix". Use -pbc mol instead.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
OUTDIR="${1:?Usage: bash process_trajectory.sh OUTDIR [REP]}"
REP="${2:-000}"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
# Unlike sbatch scripts (which SLURM copies to a temp path), regular scripts
# can reliably use BASH_SOURCE[0] to find their own location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  set +u; source "$GMXRC"; set -u
fi

if command -v gmx_mpi &>/dev/null; then
  GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then
  GMX="gmx"
else
  echo "[ERROR] No GROMACS binary found. Source your GMXRC or set GROMACS_SCRIPTS_DIR."
  exit 1
fi

# ── Input files ───────────────────────────────────────────────────────────────
TPR="${OUTDIR}/prod/rep${REP}/remd.tpr"
XTC="${OUTDIR}/trajectories/remd_rep${REP}.xtc"

[[ -f "$TPR" ]] || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -f "$XTC" ]] || { echo "[ERROR] XTC not found: $XTC"; exit 1; }

# ── Output ────────────────────────────────────────────────────────────────────
ANALYSIS_DIR="${OUTDIR}/analysis"
mkdir -p "$ANALYSIS_DIR"
OUT_XTC="${ANALYSIS_DIR}/remd_rep${REP}_pbc.xtc"

echo "[INFO] Input TPR: $TPR"
echo "[INFO] Input XTC: $XTC"
echo "[INFO] Output:    $OUT_XTC"
echo ""

# ── PBC correction ────────────────────────────────────────────────────────────
# Select "Protein" for centering, "System" for output.
printf "Protein\nSystem\n" | $GMX trjconv \
  -s "$TPR" \
  -f "$XTC" \
  -o "$OUT_XTC" \
  -pbc mol \
  -center \
  -ur compact

echo ""
echo "[OK] PBC-corrected trajectory written to: $OUT_XTC"
echo ""
echo "Next steps:"
echo "  Strip waters/ions:  bash strip_trajectory.sh $OUTDIR $REP"
echo "  Align to reference: use gmx trjconv -fit rot+trans"
