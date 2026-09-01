#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# multichain_extract_solvated_snapshots.sh — solvated snapshots for a COMPLEX
# ─────────────────────────────────────────────────────────────────────────────
# Multi-chain counterpart of extract_solvated_snapshots.sh. Writes a handful of
# PBC-corrected PDBs holding the complex plus every whole water/ion within
# SHELL_NM of it, all fitted onto a common reference so they superimpose.
#
# The ONLY difference from the single-chain script is the PBC step. Plain
# `-pbc mol` wraps each chain's centre of mass into the box *independently*, so the
# chains can land in different periodic images and the complex is split across a
# boundary — the shell would then be drawn around a torn interface, which looks
# plausible and is wrong. Chains are kept in ONE image with three per-frame passes
# (same sequence as multichain_fix_PBC.sh):
#   1. -pbc whole    : make every molecule whole.
#   2. -pbc cluster  : pull the protein chains into the same periodic image
#                      (cluster group = Protein). NOTE: "-pbc cluster" is a
#                      PERIODIC-IMAGE operation — it is NOT conformational
#                      clustering (gmx cluster); this script handles PBC only.
#   3. -pbc mol -center -ur compact : centre the complex, compact box for viewing.
# All three are per-frame (no frame-to-frame comparison).
#
# ASSUMES the complex stays within ~half the (minimum) box vector. If chains
# dissociate further the periodic image is genuinely ambiguous (and the box is too
# small — a minimum-image violation). run_analysis emits an inter-chain minimum-
# distance curve so that case is visible.
#
# Order of operations is PBC → select → fit → write, for the reasons documented in
# extract_solvated_snapshots.sh: why these are standalone PDBs rather than a
# trajectory, why the PBC passes must precede the shell selection, and why the
# selection must be computed BEFORE -fit (fitting rotates coordinates but not the
# box, so a PBC-aware `within` on a fitted frame picks up spurious far-away waters).
#
# Usage:
#   bash multichain_extract_solvated_snapshots.sh TPR XTC OUT_PREFIX [SHELL_NM] [N_SNAPSHOTS]
#
#   TPR           run input matching the XTC atom count (full system)
#   XTC           raw trajectory
#   OUT_PREFIX    output path stem; writes <OUT_PREFIX>_solvshell_<NNNNNN>ps.pdb
#   SHELL_NM      solvent shell thickness in nm (default 0.5 = 5 Å — GROMACS is nm)
#   N_SNAPSHOTS   how many evenly-spaced frames (default 5: first, 3 middle, last)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
USAGE="Usage: bash multichain_extract_solvated_snapshots.sh TPR XTC OUT_PREFIX [SHELL_NM] [N_SNAPSHOTS]"
TPR="${1:?$USAGE}"
XTC="${2:?$USAGE}"
OUT_PREFIX="${3:?$USAGE}"
# NOTE: SHELL_NM, never SHELL — $SHELL is the user's login shell.
SHELL_NM="${4:-0.5}"
N_SNAPSHOTS="${5:-5}"

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

# The loop rewrites the same temp paths once per snapshot; without this GROMACS
# squirrels away a '#md.tmp_pbc.gro.1#' backup of each, which the EXIT trap does not
# match and which would litter the output dir. Every output here is regenerable.
export GMX_MAXBACKUP=-1

# ── Validate inputs ───────────────────────────────────────────────────────────
[[ -f "$TPR" ]] || { echo "[ERROR] TPR not found: $TPR"; exit 1; }
[[ -e "$XTC" ]] || { echo "[ERROR] XTC not found: $XTC"; exit 1; }
(( N_SNAPSHOTS >= 2 )) || { echo "[ERROR] N_SNAPSHOTS must be >= 2 (got: $N_SNAPSHOTS)"; exit 1; }
awk -v s="$SHELL_NM" 'BEGIN{ exit !(s+0 > 0) }' \
  || { echo "[ERROR] SHELL_NM must be > 0 nm (got: $SHELL_NM)"; exit 1; }

mkdir -p "$(dirname "$OUT_PREFIX")"

echo "[INFO] Input TPR:   $TPR"
echo "[INFO] Input XTC:   $XTC"
echo "[INFO] Out prefix:  $OUT_PREFIX"
echo "[INFO] Shell:       ${SHELL_NM} nm"
echo "[INFO] Snapshots:   ${N_SNAPSHOTS}"
echo ""

# ── Temporaries (single frames; removed even on error) ────────────────────────
TMP_LAST="${OUT_PREFIX}.tmp_lastframe.gro"
TMP_W="${OUT_PREFIX}.tmp_whole.gro"
TMP_C="${OUT_PREFIX}.tmp_clust.gro"
TMP_PBC="${OUT_PREFIX}.tmp_pbc.gro"
TMP_FIT="${OUT_PREFIX}.tmp_fit.gro"
TMP_NDX="${OUT_PREFIX}.tmp_shell.ndx"
REF_FULL="${OUT_PREFIX}.tmp_ref_fullsys.gro"
trap 'rm -f "$TMP_LAST" "$TMP_W" "$TMP_C" "$TMP_PBC" "$TMP_FIT" "$TMP_NDX" "$REF_FULL"' EXIT

# ── Trajectory time range ─────────────────────────────────────────────────────
# Dump the final frame and read the time GROMACS stamps into the .gro title
# ("<title> t= 50000.00000 step= ..."). That string is written unconditionally, so it
# is the reliable source. `gmx check`'s "Last frame N time T" progress readout is NOT —
# it is throttled and on a long trajectory is never printed at all (see the note in
# extract_solvated_snapshots.sh), which makes the parse fail for only *some* inputs.
echo "[CMD] printf 'System\\n' | $GMX trjconv -s $TPR -f $XTC -o $TMP_LAST -dump 999999999"
printf "System\n" | $GMX trjconv \
  -s "$TPR" -f "$XTC" -o "$TMP_LAST" -dump 999999999

LAST_PS=$(head -1 "$TMP_LAST" | sed -n 's/.*t=[[:space:]]*\([0-9.eE+-]\{1,\}\).*/\1/p')
[[ -n "$LAST_PS" ]] || {
  echo "[ERROR] Could not read the final frame time from the .gro title of: $TMP_LAST"
  echo "        title was: $(head -1 "$TMP_LAST")"; exit 1; }
echo "[INFO] Last frame:  ${LAST_PS} ps"

mapfile -t TIMES < <(awk -v L="$LAST_PS" -v n="$N_SNAPSHOTS" \
  'BEGIN{ for (i=0; i<n; i++) printf "%.3f\n", L*i/(n-1) }')
echo "[INFO] Times (ps):  ${TIMES[*]}"
echo ""

# The fit reference must be FULL SYSTEM: gmx trjconv -s must match the frame's atom
# count, or the Jacobi rotation fit fails (a protein-only ref silently truncates).

for t in "${TIMES[@]}"; do
  LABEL=$(awk -v t="$t" 'BEGIN{ printf "%06.0f", t }')
  OUT_PDB="${OUT_PREFIX}_solvshell_${LABEL}ps.pdb"
  echo "─── t = ${t} ps → $(basename "$OUT_PDB")"

  # 1a. Make every molecule whole.
  echo "[CMD] printf 'System\\n' | $GMX trjconv -s $TPR -f $XTC -o $TMP_W -dump $t -pbc whole"
  printf "System\n" | $GMX trjconv \
    -s "$TPR" -f "$XTC" -o "$TMP_W" -dump "$t" -pbc whole

  # 1b. Bring the protein chains into one periodic image.
  echo "[CMD] printf 'Protein\\nSystem\\n' | $GMX trjconv -s $TPR -f $TMP_W -o $TMP_C -pbc cluster"
  printf "Protein\nSystem\n" | $GMX trjconv \
    -s "$TPR" -f "$TMP_W" -o "$TMP_C" -pbc cluster

  # 1c. Centre the complex, compact cell.
  echo "[CMD] printf 'Protein\\nSystem\\n' | $GMX trjconv -s $TPR -f $TMP_C -o $TMP_PBC -pbc mol -center -ur compact"
  printf "Protein\nSystem\n" | $GMX trjconv \
    -s "$TPR" -f "$TMP_C" -o "$TMP_PBC" -pbc mol -center -ur compact

  # First snapshot doubles as the common alignment reference.
  [[ -f "$REF_FULL" ]] || cp "$TMP_PBC" "$REF_FULL"

  # 2. Dynamic shell → index, computed on the UNFITTED frame (see the header note:
  #    -fit rotates coordinates but leaves the box vectors alone, so a PBC-aware
  #    `within` on a fitted frame measures against a box that no longer matches).
  #    `same residue as` keeps WHOLE waters/ions; without it the selection would cut
  #    molecules in half at the shell boundary. A leading quoted string names a
  #    selection; `Shell = ...` would instead declare a *variable* and gmx select
  #    would exit with "Too few selections provided".
  echo "[CMD] $GMX select -s $TPR -f $TMP_PBC -on $TMP_NDX -select '\"Shell\" group \"Protein\" or same residue as (within ${SHELL_NM} of group \"Protein\")'"
  $GMX select \
    -s "$TPR" -f "$TMP_PBC" -on "$TMP_NDX" \
    -select "\"Shell\" group \"Protein\" or same residue as (within ${SHELL_NM} of group \"Protein\")"

  # 3. Fit onto the common reference so all snapshots superimpose. Output group
  #    System, so the solvent rides along through the same rotation/translation.
  #    Atom numbering is unchanged, so the index from step 2 still applies.
  echo "[CMD] printf 'Backbone\\nSystem\\n' | $GMX trjconv -s $REF_FULL -f $TMP_PBC -o $TMP_FIT -fit rot+trans"
  printf "Backbone\nSystem\n" | $GMX trjconv \
    -s "$REF_FULL" -f "$TMP_PBC" -o "$TMP_FIT" \
    -fit rot+trans

  # 4. Write the snapshot. Select group by INDEX 0, not by name: for a dynamic
  #    selection gmx select stamps the frame/time into the group name (it becomes
  #    e.g. "Shell_f0_t1000.000"), so a name match would fail. One selection over one
  #    frame ⇒ the ndx holds exactly one group, so 0 is unambiguous.
  echo "[CMD] printf '0\\n' | $GMX trjconv -s $TPR -f $TMP_FIT -n $TMP_NDX -o $OUT_PDB"
  printf "0\n" | $GMX trjconv \
    -s "$TPR" -f "$TMP_FIT" -n "$TMP_NDX" -o "$OUT_PDB"

  [[ -s "$OUT_PDB" ]] || { echo "[ERROR] Snapshot is empty: $OUT_PDB"; exit 1; }
  echo ""
done

echo "[OK] ${N_SNAPSHOTS} solvated snapshots (${SHELL_NM} nm shell, chains kept together) written to: $(dirname "$OUT_PREFIX")/"
