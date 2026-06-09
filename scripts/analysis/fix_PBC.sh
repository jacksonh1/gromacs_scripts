#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# fix_PBC.sh — fix PBC artifacts in a trajectory
# ─────────────────────────────────────────────────────────────────────────────
# Layout-blind: takes explicit input/output paths so the same script works for
# any pipeline (plain MD, T-REMD, …). The caller supplies the concrete paths.
#
# Usage:
#   bash fix_PBC.sh TPR XTC OUT_PBC_XTC
#
#   TPR          run input matching the XTC atom count (full system)
#   XTC          raw trajectory to correct
#   OUT_PBC_XTC  output path for the PBC-corrected (full-system) trajectory
#
# PBC treatment: -pbc mol -center -ur compact
#   - mol:     makes each molecule whole and puts its COM inside the box.
#              Per-frame operation — safe for REMD trajectories.
#   - center:  centers the protein in the box.
#   - compact: renders the dodecahedron box as a compact shape; without this
#              the box looks "exploded" in VMD/PyMOL.
#
# This is the SINGLE-CHAIN path. For a multi-chain complex, -pbc mol wraps each
# chain's COM independently and can split the complex across a box boundary — use
# the multichain_* analysis scripts instead (run_analysis.sh dispatches on chain
# count automatically).
#
# WARNING: -pbc nojump is intentionally NOT used here. nojump compares
# consecutive frames to detect box-boundary crossings, but in T-REMD,
# coordinate exchanges cause discontinuous jumps between frames that nojump
# would misinterpret and incorrectly "fix". Use -pbc mol instead.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
TPR="${1:?Usage: bash fix_PBC.sh TPR XTC OUT_PBC_XTC}"
XTC="${2:?Usage: bash fix_PBC.sh TPR XTC OUT_PBC_XTC}"
OUT_XTC="${3:?Usage: bash fix_PBC.sh TPR XTC OUT_PBC_XTC}"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
# Unlike sbatch scripts (which SLURM copies to a temp path), regular scripts
# can reliably use BASH_SOURCE[0] to find their own location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  set +u; source "$GMXRC"; set -u
fi

if command -v gmx_mpi &>/dev/null; then
  GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then
  GMX="gmx"
else
  echo "[ERROR] No GROMACS binary (gmx_mpi/gmx) on PATH. Source your GROMACS GMXRC or load the GROMACS module first."
  exit 1
fi

# ── Validate inputs ───────────────────────────────────────────────────────────
[[ -f "$TPR" ]] || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -f "$XTC" ]] || { echo "[ERROR] XTC not found: $XTC"; exit 1; }

mkdir -p "$(dirname "$OUT_XTC")"

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
