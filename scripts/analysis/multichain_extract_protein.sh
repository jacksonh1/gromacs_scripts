#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_extract_protein.sh — protein-only RMSD reference for a complex
# ─────────────────────────────────────────────────────────────────────────────
# Multi-chain counterpart of extract_protein.sh. Selects the "Protein" group from
# a full-system structure (e.g. em/em.gro) and writes it on its own, using
# `-pbc cluster` so the chains are whole AND kept in the SAME periodic image.
#
# Why cluster and not whole: `gmx rms` does no PBC treatment of its reference, so a
# complex left split across a box boundary (one chain a box-vector away) silently
# corrupts every RMSD value. `-pbc whole` only un-breaks within a molecule; it does
# not bring separate chains together. NOTE: "-pbc cluster" is a periodic-image
# operation, NOT conformational clustering (gmx cluster).
#
# Output atom ordering matches the protein-only trajectory from
# strip_and_align_trajectory.sh, so it works directly as the reference for the
# whole-complex and per-chain RMSD.
#
# Usage:
#   bash multichain_extract_protein.sh STRUCT TPR OUT_GRO
#
#   STRUCT   full-system structure (.gro/.pdb), e.g. em/em.gro
#   TPR      run input matching STRUCT's atom count (provides the "Protein" group)
#   OUT_GRO  protein-only output (.gro/.pdb)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

STRUCT="${1:?Usage: bash multichain_extract_protein.sh STRUCT TPR OUT_GRO}"
TPR="${2:?Usage: bash multichain_extract_protein.sh STRUCT TPR OUT_GRO}"
OUT_GRO="${3:?Usage: bash multichain_extract_protein.sh STRUCT TPR OUT_GRO}"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_CONFIG="${SCRIPT_DIR}/../../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  set +u; source "$GMXRC"; set -u
fi

if command -v gmx_mpi &>/dev/null; then GMX="gmx_mpi"
elif command -v gmx &>/dev/null; then GMX="gmx"
else echo "[ERROR] No GROMACS binary (gmx_mpi/gmx) on PATH. Source your GROMACS GMXRC or load the GROMACS module first."; exit 1
fi

[[ -f "$STRUCT" ]] || { echo "[ERROR] Structure not found: $STRUCT"; exit 1; }
[[ -f "$TPR" ]]    || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
mkdir -p "$(dirname "$OUT_GRO")"

echo "[INFO] Extracting Protein group (chains kept together) from $STRUCT → $OUT_GRO"
# cluster prompts twice: (1) group to cluster, (2) group to output — both Protein.
printf "Protein\nProtein\n" | $GMX trjconv -s "$TPR" -f "$STRUCT" -pbc cluster -o "$OUT_GRO"

echo "[OK] Protein-only reference written to: $OUT_GRO"
