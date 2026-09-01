#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# finalize_outputs.sh — make the scratch run directory a self-contained record
# ─────────────────────────────────────────────────────────────────────────────
# Under the folder-symlink output model, the bulk stage dirs (prod/ equil/
# density/ heat/ relax/) are OUTDIR *symlinks into scratch*, so they are already
# real on scratch after the run. The small, laptop-worthy dirs stay REAL in
# OUTDIR (analysis/ solvated_snapshots/ logs/ em/ build/ parameters.txt,
# input_structure.pdb, the exported final PDB). This copies those into SCRATCH_DIR at
# the same relative path so the scratch archive stands alone (self-describing:
# every scratch dir then carries its own parameters.txt + analysis/).
#
# Usage:
#   bash finalize_outputs.sh OUTDIR SCRATCH_DIR
#
# Idempotent: each destination is removed and re-copied, so re-running (or a
# restart) refreshes the archive cleanly. The stage symlinks are skipped on
# purpose — copying one would recurse back into SCRATCH_DIR.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

OUTDIR="${1:?Usage: finalize_outputs.sh OUTDIR SCRATCH_DIR}"
SCRATCH_DIR="${2:?Usage: finalize_outputs.sh OUTDIR SCRATCH_DIR}"
OUTDIR="${OUTDIR%/}"
SCRATCH_DIR="${SCRATCH_DIR%/}"

[[ -d "$OUTDIR" ]]      || { echo "[ERROR] OUTDIR not found: $OUTDIR"; exit 1; }
[[ -d "$SCRATCH_DIR" ]] || { echo "[ERROR] SCRATCH_DIR not found: $SCRATCH_DIR"; exit 1; }

# The real (non-symlink) OUTDIR entries worth archiving. topol/ is only present
# for REST2; the loop skips anything that does not exist. A stage dir that is a
# symlink (prod/ equil/ density/ heat/ relax/) is never copied — it already lives
# on scratch and copying it would recurse. Final PDB(s) are handled below.
for entry in analysis solvated_snapshots logs em build topol parameters.txt input_structure.pdb; do
  src="${OUTDIR}/${entry}"
  [[ -e "$src" ]] || continue      # not present for this engine → skip
  [[ -L "$src" ]] && continue      # a stage symlink → already on scratch, never copy
  dst="${SCRATCH_DIR}/${entry}"
  echo "[CMD] rm -rf ${dst} && cp -a ${src} ${dst}"
  rm -rf "$dst"
  cp -a "$src" "$dst"
done

# Exported final PDB(s): <base>_final.pdb (MD) or <base>_final_rep000.pdb (REMD/REST2).
shopt -s nullglob
for pdb in "${OUTDIR}"/*_final*.pdb; do
  dst="${SCRATCH_DIR}/$(basename "$pdb")"
  echo "[CMD] cp -a ${pdb} ${dst}"
  rm -rf "$dst"
  cp -a "$pdb" "$dst"
done
shopt -u nullglob

echo "[OK] Scratch archive is self-contained → $SCRATCH_DIR"
