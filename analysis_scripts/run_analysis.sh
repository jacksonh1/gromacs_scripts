#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_analysis.sh — re-run the full post-analysis for an MD or T-REMD job
# ─────────────────────────────────────────────────────────────────────────────
# One entry point for the whole analysis step, so it can be re-run by hand
# (e.g. after fixing an analysis script) without resubmitting the simulation.
# Both engines' sbatch scripts call this too, so there is a single copy of the
# orchestration to maintain.
#
# Usage:
#   bash run_analysis.sh OUTDIR [REP]
#
#   OUTDIR   a finished job directory (the one holding build/ em/ prod/ ...)
#   REP      T-REMD replica slot to analyse (default 000; ignored for plain MD)
#
# Pipeline (MD or REMD):
#   [REMD only] remd_acceptance.py                 — exchange acceptance rates
#   PBC fix + strip + align → <prefix>_stripped_aligned.{xtc,gro}
#   protein reference       → <prefix>_init.gro    — minimized RMSD reference
#   calc_traj_rmsd/rg/rmsf/dssp + plot_xvg/plot_dssp   (whole protein/complex)
#   cluster_traj.py                                — conformational clustering (Cα)
#   [multi-chain only] per-chain RMSD/RMSF + inter-chain min distance
#
# Two auto-detections, so one command serves every job:
#   - MD vs REMD     — from the job layout (prod/md.tpr vs prod/rep<REP>/remd.tpr).
#   - chain count    — from the topology. 1 chain → the simple path (fix_PBC.sh /
#                      extract_protein.sh). >1 chain → the multichain_* path, which
#                      keeps the complex's chains in one periodic image (-pbc cluster)
#                      and adds the per-chain / inter-chain metrics. Same core output
#                      files either way.
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail   # not -e: one failed metric should not abort the rest

OUTDIR="${1:?Usage: bash run_analysis.sh OUTDIR [REP]}"
REP="${2:-000}"
OUTDIR="${OUTDIR%/}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Overwrite analysis outputs cleanly on re-run instead of leaving GROMACS '#file.N#'
# backups. Safe here: every output is regenerable from the trajectory + tpr, which
# this script never writes. (Only affects the gmx calls launched below.)
export GMX_MAXBACKUP=-1

# Activate the analysis/plotting Python env (matplotlib; needed by the plotters
# and remd_acceptance.py). The per-metric calc scripts source GROMACS themselves.
SITE_CONFIG="${SCRIPT_DIR}/../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  if declare -f activate_python_env >/dev/null; then activate_python_env; fi
fi

[[ -d "$OUTDIR" ]] || { echo "[ERROR] Not a directory: $OUTDIR"; exit 1; }

# ── Auto-detect pipeline + resolve paths ──────────────────────────────────────
ANALYSIS_DIR="${OUTDIR}/analysis"
EM_GRO="${OUTDIR}/em/em.gro"
EM_TPR="${OUTDIR}/em/em.tpr"

if [[ -f "${OUTDIR}/prod/md.tpr" ]]; then
  MODE="MD"
  TPR="${OUTDIR}/prod/md.tpr"
  XTC="${OUTDIR}/trajectories/md.xtc"
  PREFIX="${ANALYSIS_DIR}/md"
elif [[ -f "${OUTDIR}/prod/rep${REP}/remd.tpr" ]]; then
  MODE="REMD"
  TPR="${OUTDIR}/prod/rep${REP}/remd.tpr"
  XTC="${OUTDIR}/trajectories/remd_rep${REP}.xtc"
  PREFIX="${ANALYSIS_DIR}/remd_rep${REP}"
else
  echo "[ERROR] Cannot find prod/md.tpr or prod/rep${REP}/remd.tpr under $OUTDIR"
  echo "        Is this a finished MD or T-REMD job directory?"
  exit 1
fi

mkdir -p "$ANALYSIS_DIR"

echo "===== run_analysis: ${MODE} ====="
echo "[INFO] OUTDIR : $OUTDIR"
echo "[INFO] TPR    : $TPR"
echo "[INFO] XTC    : $XTC"
echo "[INFO] PREFIX : $PREFIX"

[[ -f "$TPR" ]]  || { echo "[ERROR] Run input not found: $TPR"; exit 1; }
# -e follows the symlink: trajectories/ are symlinks into scratch, which is purged.
[[ -e "$XTC" ]]  || { echo "[ERROR] Trajectory not found (scratch purged?): $XTC"; exit 1; }

# ── REMD-only: exchange acceptance rates ──────────────────────────────────────
if [[ "$MODE" == "REMD" ]]; then
  echo "[CMD] python3 ${SCRIPT_DIR}/remd_acceptance.py $OUTDIR"
  python3 "${SCRIPT_DIR}/remd_acceptance.py" "$OUTDIR" \
    || echo "[WARN] remd_acceptance.py failed — re-run the command above"
fi

# ── Detect protein chain count (topology) → single- vs multi-chain path ───────
# Count protein molecules in the system topology's [ molecules ] section. >1 chain
# ⇒ the complex must be kept in one periodic image (the multichain_* pipeline).
BUILD_DIR="${OUTDIR}/build"
NCHAINS=1
TOP=$(ls "${BUILD_DIR}"/*.top 2>/dev/null | head -1 || true)
if [[ -n "$TOP" ]]; then
  NCHAINS=$(awk '
    /^[[:space:]]*\[/ { s=$0; gsub(/[][[:space:]]/,"",s); insec=(tolower(s)=="molecules"); next }
    insec { line=$0; sub(/;.*/,"",line); n=split(line,a," ");
            if (n>=2 && a[1] ~ /^Protein/) tot+=a[2] }
    END { print (tot>0 ? tot : 1) }
  ' "$TOP")
fi

# ── 1. PBC fix + strip + align (+ RMSD reference) → <prefix>_stripped_aligned.* ─
# Dispatch on chain count: single chain keeps the simple path; a multi-chain
# complex uses the multichain_* scripts (which keep the chains in one image). Both
# produce the SAME output files, so everything downstream is identical.
INIT_REF="${PREFIX}_init.gro"
if (( NCHAINS > 1 )); then
  echo "[INFO] Multi-chain system (${NCHAINS} protein chains) → multichain pipeline"
  echo "[CMD] bash ${SCRIPT_DIR}/multichain_fix_PBC_strip_align.sh $TPR $XTC $PREFIX"
  bash "${SCRIPT_DIR}/multichain_fix_PBC_strip_align.sh" "$TPR" "$XTC" "$PREFIX" \
    || echo "[WARN] multichain_fix_PBC_strip_align.sh failed — re-run the command above"
  echo "[CMD] bash ${SCRIPT_DIR}/multichain_extract_protein.sh $EM_GRO $EM_TPR $INIT_REF"
  bash "${SCRIPT_DIR}/multichain_extract_protein.sh" "$EM_GRO" "$EM_TPR" "$INIT_REF" \
    || echo "[WARN] multichain_extract_protein.sh failed — RMSD will fall back to frame-0 reference"
else
  echo "[INFO] Single-chain system → standard pipeline"
  echo "[CMD] bash ${SCRIPT_DIR}/fix_PBC_strip_align.sh $TPR $XTC $PREFIX"
  bash "${SCRIPT_DIR}/fix_PBC_strip_align.sh" "$TPR" "$XTC" "$PREFIX" \
    || echo "[WARN] fix_PBC_strip_align.sh failed — re-run the command above"
  echo "[CMD] bash ${SCRIPT_DIR}/extract_protein.sh $EM_GRO $EM_TPR $INIT_REF"
  bash "${SCRIPT_DIR}/extract_protein.sh" "$EM_GRO" "$EM_TPR" "$INIT_REF" \
    || echo "[WARN] extract_protein.sh failed — RMSD will fall back to frame-0 reference"
fi

REF="${PREFIX}_stripped_aligned.gro"   # first frame; topology for Rg/RMSF/DSSP
FIT="${PREFIX}_stripped_aligned.xtc"
# RMSD reference is the minimized starting structure (drift from the design); pass
# "$REF" instead for frame-0 RMSD. Fall back to frame-0 if extraction failed.
[[ -f "$INIT_REF" ]] || INIT_REF="$REF"

# ── 2. Metrics (on the stripped/aligned protein trajectory) + plots ───────────
echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_rmsd.sh $INIT_REF $FIT ${PREFIX}_rmsd.xvg"
bash "${SCRIPT_DIR}/calc_traj_rmsd.sh" "$INIT_REF" "$FIT" "${PREFIX}_rmsd.xvg" \
  && python3 "${SCRIPT_DIR}/plot_xvg.py" "${PREFIX}_rmsd.xvg" "${PREFIX}_rmsd.png" \
  || echo "[WARN] RMSD step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_rg.sh $REF $FIT ${PREFIX}_rg.xvg"
bash "${SCRIPT_DIR}/calc_traj_rg.sh" "$REF" "$FIT" "${PREFIX}_rg.xvg" \
  && python3 "${SCRIPT_DIR}/plot_xvg.py" "${PREFIX}_rg.xvg" "${PREFIX}_rg.png" \
  || echo "[WARN] Rg step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_rmsf.sh $REF $FIT ${PREFIX}_rmsf.xvg"
bash "${SCRIPT_DIR}/calc_traj_rmsf.sh" "$REF" "$FIT" "${PREFIX}_rmsf.xvg" \
  && python3 "${SCRIPT_DIR}/plot_xvg.py" "${PREFIX}_rmsf.xvg" "${PREFIX}_rmsf.png" \
  || echo "[WARN] RMSF step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_dssp.sh $REF $FIT ${PREFIX}_dssp.dat"
bash "${SCRIPT_DIR}/calc_traj_dssp.sh" "$REF" "$FIT" "${PREFIX}_dssp.dat" \
  && python3 "${SCRIPT_DIR}/plot_dssp.py" "${PREFIX}_dssp.dat" "${PREFIX}_dssp.png" \
  || echo "[WARN] DSSP step failed — re-run the commands above"

# Conformational clustering (Cα-RMSD; sklearn DBSCAN). Runs on the same protein-only
# aligned trajectory, so it serves single- and multi-chain alike. CLUSTER_CUTOFF is
# the backbone-RMSD cutoff in nm (default 0.20 = 2.0 Å); see cluster_traj.py.
echo "[CMD] python3 ${SCRIPT_DIR}/cluster_traj.py $REF $FIT ${PREFIX} --cutoff ${CLUSTER_CUTOFF:-0.20}"
python3 "${SCRIPT_DIR}/cluster_traj.py" "$REF" "$FIT" "${PREFIX}" --cutoff "${CLUSTER_CUTOFF:-0.20}" \
  || echo "[WARN] clustering step failed — re-run the command above"

# ── 3. Multi-chain extras: per-chain RMSD/RMSF + inter-chain min distance ──────
if (( NCHAINS > 1 )); then
  echo ""
  echo "[INFO] Multi-chain extras (per-chain RMSD/RMSF + inter-chain distance)"
  NDX="${PREFIX}_chains.ndx"
  echo "[CMD] python3 ${SCRIPT_DIR}/multichain_chain_index.py $BUILD_DIR $REF $NDX"
  IDX_OUT=$(python3 "${SCRIPT_DIR}/multichain_chain_index.py" "$BUILD_DIR" "$REF" "$NDX" 2>&1) \
    || echo "[WARN] multichain_chain_index.py failed"
  echo "$IDX_OUT"
  CHAINS=( $(printf '%s\n' "$IDX_OUT" | sed -n 's/^CHAINS: //p') )

  if (( ${#CHAINS[@]} >= 2 )); then
    # Per-chain backbone RMSD (vs init reference) and per-residue RMSF.
    for c in "${CHAINS[@]}"; do
      g="Chain${c}_Backbone"
      echo "[CMD] bash ${SCRIPT_DIR}/multichain_chain_rmsd.sh $INIT_REF $FIT $NDX $g ${PREFIX}_chain${c}_rmsd.xvg"
      bash "${SCRIPT_DIR}/multichain_chain_rmsd.sh" "$INIT_REF" "$FIT" "$NDX" "$g" "${PREFIX}_chain${c}_rmsd.xvg" \
        && python3 "${SCRIPT_DIR}/plot_xvg.py" "${PREFIX}_chain${c}_rmsd.xvg" "${PREFIX}_chain${c}_rmsd.png" \
        || echo "[WARN] chain ${c} RMSD failed"
      echo "[CMD] bash ${SCRIPT_DIR}/multichain_chain_rmsf.sh $REF $FIT $NDX $g ${PREFIX}_chain${c}_rmsf.xvg"
      bash "${SCRIPT_DIR}/multichain_chain_rmsf.sh" "$REF" "$FIT" "$NDX" "$g" "${PREFIX}_chain${c}_rmsf.xvg" \
        && python3 "${SCRIPT_DIR}/plot_xvg.py" "${PREFIX}_chain${c}_rmsf.xvg" "${PREFIX}_chain${c}_rmsf.png" \
        || echo "[WARN] chain ${c} RMSF failed"
    done
    # Inter-chain minimum distance for each unique chain pair (binding observable).
    n=${#CHAINS[@]}
    for ((i=0; i<n; i++)); do
      for ((j=i+1; j<n; j++)); do
        a="${CHAINS[i]}"; b="${CHAINS[j]}"
        if (( n == 2 )); then out="${PREFIX}_interchain_mindist"; else out="${PREFIX}_interchain_${a}_${b}_mindist"; fi
        echo "[CMD] bash ${SCRIPT_DIR}/multichain_interchain_dist.sh $REF $FIT $NDX Chain${a} Chain${b} ${out}.xvg"
        bash "${SCRIPT_DIR}/multichain_interchain_dist.sh" "$REF" "$FIT" "$NDX" "Chain${a}" "Chain${b}" "${out}.xvg" \
          && python3 "${SCRIPT_DIR}/plot_xvg.py" "${out}.xvg" "${out}.png" \
          || echo "[WARN] inter-chain ${a}-${b} distance failed"
      done
    done
  else
    echo "[WARN] chain index produced <2 chains; skipping per-chain extras"
  fi
fi

echo "[OK] Post-analysis done → ${ANALYSIS_DIR}/"
