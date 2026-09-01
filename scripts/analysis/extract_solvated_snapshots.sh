#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# extract_solvated_snapshots.sh — a few PBC-corrected PDBs WITH a solvent shell
# ─────────────────────────────────────────────────────────────────────────────
# Every other analysis output is protein-only (strip_and_align_trajectory.sh drops
# waters/ions to keep the aligned trajectory small). That is right for RMSD/Rg/RMSF,
# but it hides the solvent that mediates a binding interface. This writes a HANDFUL
# of snapshots that keep the protein plus every whole water/ion within SHELL_NM of
# it — evenly spaced across the trajectory, all fitted onto a common reference so
# they superimpose when loaded together.
#
# Layout-blind: explicit input paths + an output prefix (same contract as fix_PBC.sh).
#
# Usage:
#   bash extract_solvated_snapshots.sh TPR XTC OUT_PREFIX [SHELL_NM] [N_SNAPSHOTS]
#
#   TPR           run input matching the XTC atom count (full system)
#   XTC           raw trajectory
#   OUT_PREFIX    output path stem; writes <OUT_PREFIX>_solvshell_<NNNNNN>ps.pdb
#   SHELL_NM      solvent shell thickness in nm (default 0.5 = 5 Å — GROMACS is nm)
#   N_SNAPSHOTS   how many evenly-spaced frames (default 5: first, 3 middle, last)
#
# WHY A HANDFUL OF PDBs AND NOT A TRAJECTORY
# A shell selection is *dynamic*: the set of waters within SHELL_NM of the protein
# changes every frame, so the atom count changes every frame. XTC/TRR require a
# fixed atom count and there is no single topology that matches such a trajectory.
# Independent one-frame PDBs sidestep that entirely — each carries its own atoms.
#
# ORDER OF OPERATIONS: PBC → select → fit → write. All three constraints are real.
#   1. PBC before selection. `gmx select`'s `within` is PBC-aware, so it will happily
#      select a water whose *periodic image* is near the protein while its stored
#      coordinate sits across the box — written into the PDB as a water floating in
#      space. -pbc mol -center -ur compact first puts the protein at the box centre
#      with all nearest images already in place.
#   2. SELECT BEFORE FITTING. -fit rot+trans rotates the coordinates but leaves the
#      box vectors untouched, so on a fitted frame the PBC-aware `within` measures
#      against a box that no longer corresponds to the coordinates and picks up
#      spurious far-away waters (measured: 1622 atoms selected on the fitted frame
#      vs 1598 on the same frame pre-fit — 8 bogus waters up to 23 Å out). Atom
#      numbering is identical either way, so the index applies to the fitted frame.
#   3. PBC and fitting must be SEPARATE trjconv calls. trjconv refuses -fit together
#      with -pbc mol (this is also why fix_PBC.sh and strip_and_align_trajectory.sh
#      are separate scripts).
#
# This is the SINGLE-CHAIN path. For a multi-chain complex -pbc mol wraps each
# chain's COM independently and can split the complex — a "shell" around a torn
# interface. Use multichain_extract_solvated_snapshots.sh instead (run_analysis.sh
# dispatches on chain count automatically).
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
USAGE="Usage: bash extract_solvated_snapshots.sh TPR XTC OUT_PREFIX [SHELL_NM] [N_SNAPSHOTS]"
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
TMP_PBC="${OUT_PREFIX}.tmp_pbc.gro"
TMP_FIT="${OUT_PREFIX}.tmp_fit.gro"
TMP_NDX="${OUT_PREFIX}.tmp_shell.ndx"
REF_FULL="${OUT_PREFIX}.tmp_ref_fullsys.gro"
trap 'rm -f "$TMP_LAST" "$TMP_PBC" "$TMP_FIT" "$TMP_NDX" "$REF_FULL"' EXIT

# ── Trajectory time range ─────────────────────────────────────────────────────
# Dump the final frame and read the time GROMACS stamps into the .gro title
# ("<title> t= 50000.00000 step= ..."). That string is written unconditionally, so it
# is the reliable source.
#
# Do NOT parse `gmx check` for this. Its "Last frame N time T" progress readout is
# THROTTLED: on a long trajectory it is never printed at all (measured on a 2501-frame,
# 808 MB xtc — progress stopped at "Reading frame 2000" and no "Last frame" line was
# emitted, while gmx check still exited 0). The parse then silently yields nothing, and
# it fails only for *some* inputs, which is the worst possible failure mode.
echo "[CMD] printf 'System\\n' | $GMX trjconv -s $TPR -f $XTC -o $TMP_LAST -dump 999999999"
printf "System\n" | $GMX trjconv \
  -s "$TPR" -f "$XTC" -o "$TMP_LAST" -dump 999999999

LAST_PS=$(head -1 "$TMP_LAST" | sed -n 's/.*t=[[:space:]]*\([0-9.eE+-]\{1,\}\).*/\1/p')
[[ -n "$LAST_PS" ]] || {
  echo "[ERROR] Could not read the final frame time from the .gro title of: $TMP_LAST"
  echo "        title was: $(head -1 "$TMP_LAST")"; exit 1; }
echo "[INFO] Last frame:  ${LAST_PS} ps"

# Evenly spaced over [0, LAST_PS]. -dump snaps to the nearest stored frame, so an
# inexact division is harmless.
mapfile -t TIMES < <(awk -v L="$LAST_PS" -v n="$N_SNAPSHOTS" \
  'BEGIN{ for (i=0; i<n; i++) printf "%.3f\n", L*i/(n-1) }')
echo "[INFO] Times (ps):  ${TIMES[*]}"
echo ""

# The fit reference must be FULL SYSTEM: gmx trjconv -s must match the frame's atom
# count, or the Jacobi rotation fit fails (a protein-only ref silently truncates).
# Frame 0 of this same trajectory, after PBC correction, is that reference.

for t in "${TIMES[@]}"; do
  LABEL=$(awk -v t="$t" 'BEGIN{ printf "%06.0f", t }')
  OUT_PDB="${OUT_PREFIX}_solvshell_${LABEL}ps.pdb"
  echo "─── t = ${t} ps → $(basename "$OUT_PDB")"

  # 1. PBC: whole molecules, protein centred, compact cell (full system out).
  echo "[CMD] printf 'Protein\\nSystem\\n' | $GMX trjconv -s $TPR -f $XTC -o $TMP_PBC -dump $t -pbc mol -center -ur compact"
  printf "Protein\nSystem\n" | $GMX trjconv \
    -s "$TPR" -f "$XTC" -o "$TMP_PBC" \
    -dump "$t" -pbc mol -center -ur compact

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

echo "[OK] ${N_SNAPSHOTS} solvated snapshots (${SHELL_NM} nm shell) written to: $(dirname "$OUT_PREFIX")/"
