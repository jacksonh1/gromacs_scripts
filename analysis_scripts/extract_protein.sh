#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# extract_protein.sh — write a protein-only structure from a full-system frame
# ─────────────────────────────────────────────────────────────────────────────
# Layout-blind helper, shared by the plain-MD and T-REMD pipelines. Selects the
# "Protein" group from a full-system structure and writes it on its own. The
# output has the same atom ordering as the protein-only trajectory produced by
# strip_and_align_trajectory.sh (both select "Protein" from the same system), so
# it can be used directly as the reference for calc_traj_rmsd.sh etc.
#
# Usage:
#   bash extract_protein.sh STRUCT TPR OUT_GRO
#
#   STRUCT   full-system structure (.gro/.pdb), e.g. em/em.gro
#   TPR      run input matching STRUCT's atom count (provides the "Protein" group)
#   OUT_GRO  protein-only output (.gro/.pdb)
#
# Typical use: build an initial-structure reference from the minimized structure
# so RMSD is measured as drift from the starting structure rather than from the
# first frame of the trajectory. The protein is made whole (-pbc whole) so a
# molecule broken across the box boundary can't corrupt the RMSD reference.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash extract_protein.sh STRUCT TPR OUT_GRO}"
TPR="${2:?Usage: bash extract_protein.sh STRUCT TPR OUT_GRO}"
OUT_GRO="${3:?Usage: bash extract_protein.sh STRUCT TPR OUT_GRO}"

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

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TPR" ]]    || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
mkdir -p "$(dirname "$OUT_GRO")"

echo "[INFO] Extracting Protein group from $STRUCT → $OUT_GRO"
# -pbc whole un-breaks the protein if it straddles a box boundary in STRUCT.
# gmx rms does NOT make its reference whole, so a broken reference would silently
# corrupt every RMSD value — make it whole here. No-op if already intact.
#
# This is the SINGLE-CHAIN path. A multi-chain complex needs -pbc cluster to keep
# the chains in one image — use multichain_extract_protein.sh (run_analysis.sh
# dispatches on chain count automatically).
printf "Protein\n" | $GMX trjconv -s "$TPR" -f "$STRUCT" -pbc whole -o "$OUT_GRO"

echo "[OK] Protein-only structure written to: $OUT_GRO"
