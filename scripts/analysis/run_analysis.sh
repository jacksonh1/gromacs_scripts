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
#   [REMD only] gromd-acceptance                   — exchange acceptance rates
#   PBC fix + strip + align → <prefix>_stripped_aligned.{xtc,gro}
#   protein reference       → <prefix>_init.gro    — minimized RMSD reference
#   calc_traj_rmsd/rg/rmsf/dssp + gromd-plot-xvg/-dssp   (whole protein/complex)
#   gromd-cluster                                  — conformational clustering (Cα)
#   [MD only] solvated snapshots → <OUTDIR>/solvated_snapshots/*.pdb
#   [multi-chain only] per-chain RMSD/RMSF + inter-chain min distance
#
# Two auto-detections, so one command serves every job. Both live in
# `gromd-layout` (gromd_analysis/layout.py), which parses OUTDIR into a JobDir:
#   - MD vs REMD vs REST2 — from the job layout (prod/md.tpr vs prod/rep<REP>/remd.tpr
#                      vs prod/rep<REP>/rest2.tpr).
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
# and gromd-acceptance). The per-metric calc scripts source GROMACS themselves.
SITE_CONFIG="${SCRIPT_DIR}/../../site_config.sh"
if [[ -f "$SITE_CONFIG" ]]; then
  source "$SITE_CONFIG"
  if declare -f activate_python_env >/dev/null; then activate_python_env; fi
fi

[[ -d "$OUTDIR" ]] || { echo "[ERROR] Not a directory: $OUTDIR"; exit 1; }

# ── Auto-detect pipeline + resolve paths ──────────────────────────────────────
# gromd-layout parses the job directory into typed, existing paths — MD vs T-REMD
# vs REST2 from the layout, protein chain count from the topology's [ molecules ]
# section, plus the legacy trajectories/ fallback — and emits them as shell
# assignments. It exits non-zero with an [ERROR] when this is not a finished job,
# so a bad OUTDIR fails here instead of halfway through the metrics.
#
# Sets: MODE TPR XTC PREFIX NCHAINS ANALYSIS_DIR EM_GRO EM_TPR BUILD_DIR
echo "[CMD] eval \"\$(gromd-layout $OUTDIR $REP)\""
LAYOUT="$(gromd-layout "$OUTDIR" "$REP")" || exit 1
eval "$LAYOUT"

mkdir -p "$ANALYSIS_DIR"

echo "===== run_analysis: ${MODE} ====="
echo "[INFO] OUTDIR : $OUTDIR"
echo "[INFO] TPR    : $TPR"
echo "[INFO] XTC    : $XTC"
echo "[INFO] PREFIX : $PREFIX"

# ── REMD/REST2-only: exchange acceptance rates ────────────────────────────────
# (gromd-acceptance auto-detects the remd.log vs rest2.log basename.)
if [[ "$MODE" == "REMD" || "$MODE" == "REST2" ]]; then
  echo "[CMD] gromd-acceptance $OUTDIR"
  gromd-acceptance "$OUTDIR" \
    || echo "[WARN] gromd-acceptance failed — re-run the command above"
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
  && gromd-plot-xvg "${PREFIX}_rmsd.xvg" "${PREFIX}_rmsd.png" \
  || echo "[WARN] RMSD step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_rg.sh $REF $FIT ${PREFIX}_rg.xvg"
bash "${SCRIPT_DIR}/calc_traj_rg.sh" "$REF" "$FIT" "${PREFIX}_rg.xvg" \
  && gromd-plot-xvg "${PREFIX}_rg.xvg" "${PREFIX}_rg.png" \
  || echo "[WARN] Rg step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_rmsf.sh $REF $FIT ${PREFIX}_rmsf.xvg"
bash "${SCRIPT_DIR}/calc_traj_rmsf.sh" "$REF" "$FIT" "${PREFIX}_rmsf.xvg" \
  && gromd-plot-xvg "${PREFIX}_rmsf.xvg" "${PREFIX}_rmsf.png" \
  || echo "[WARN] RMSF step failed — re-run the commands above"

echo "[CMD] bash ${SCRIPT_DIR}/calc_traj_dssp.sh $REF $FIT ${PREFIX}_dssp.dat"
bash "${SCRIPT_DIR}/calc_traj_dssp.sh" "$REF" "$FIT" "${PREFIX}_dssp.dat" \
  && gromd-plot-dssp "${PREFIX}_dssp.dat" "${PREFIX}_dssp.png" \
  || echo "[WARN] DSSP step failed — re-run the commands above"

# Conformational clustering (Cα-RMSD; sklearn DBSCAN). Runs on the same protein-only
# aligned trajectory, so it serves single- and multi-chain alike. CLUSTER_CUTOFF is
# the backbone-RMSD cutoff in nm (default 0.20 = 2.0 Å); see gromd_analysis/clustering.py.
echo "[CMD] gromd-cluster $REF $FIT ${PREFIX} --cutoff ${CLUSTER_CUTOFF:-0.20}"
gromd-cluster "$REF" "$FIT" "${PREFIX}" --cutoff "${CLUSTER_CUTOFF:-0.20}" \
  || echo "[WARN] clustering step failed — re-run the command above"

# ── 3. Solvated snapshots (MD only) → <OUTDIR>/solvated_snapshots/ ────────────
# A few PBC-corrected, mutually-aligned PDBs that KEEP a solvent shell, for looking
# at interface waters (the metrics above all run on the protein-only trajectory).
# MD only for now — this serves bound-state sampling, not the T-REMD objectives.
# Same chain-count dispatch as step 1. SHELL_NM/N_SNAPSHOTS are analysis knobs, so
# they are env overrides here rather than job parameters (cf. CLUSTER_CUTOFF above).
if [[ "$MODE" == "MD" ]]; then
  echo ""
  SNAP_PREFIX="${OUTDIR}/solvated_snapshots/md"
  if (( NCHAINS > 1 )); then SNAP_SCRIPT="multichain_extract_solvated_snapshots.sh"
  else                       SNAP_SCRIPT="extract_solvated_snapshots.sh"; fi
  echo "[CMD] bash ${SCRIPT_DIR}/${SNAP_SCRIPT} $TPR $XTC $SNAP_PREFIX ${SHELL_NM:-0.5} ${N_SNAPSHOTS:-5}"
  bash "${SCRIPT_DIR}/${SNAP_SCRIPT}" "$TPR" "$XTC" "$SNAP_PREFIX" "${SHELL_NM:-0.5}" "${N_SNAPSHOTS:-5}" \
    || echo "[WARN] solvated snapshots failed — re-run the command above"
fi

# ── 4. Multi-chain extras: per-chain RMSD/RMSF + inter-chain min distance ──────
if (( NCHAINS > 1 )); then
  echo ""
  echo "[INFO] Multi-chain extras (per-chain RMSD/RMSF + inter-chain distance)"
  NDX="${PREFIX}_chains.ndx"
  echo "[CMD] gromd-chain-index $BUILD_DIR $REF $NDX"
  IDX_OUT=$(gromd-chain-index "$BUILD_DIR" "$REF" "$NDX" 2>&1) \
    || echo "[WARN] gromd-chain-index failed"
  echo "$IDX_OUT"
  CHAINS=( $(printf '%s\n' "$IDX_OUT" | sed -n 's/^CHAINS: //p') )

  if (( ${#CHAINS[@]} >= 2 )); then
    # Per-chain backbone RMSD (vs init reference) and per-residue RMSF.
    for c in "${CHAINS[@]}"; do
      g="Chain${c}_Backbone"
      echo "[CMD] bash ${SCRIPT_DIR}/multichain_chain_rmsd.sh $INIT_REF $FIT $NDX $g ${PREFIX}_chain${c}_rmsd.xvg"
      bash "${SCRIPT_DIR}/multichain_chain_rmsd.sh" "$INIT_REF" "$FIT" "$NDX" "$g" "${PREFIX}_chain${c}_rmsd.xvg" \
        && gromd-plot-xvg "${PREFIX}_chain${c}_rmsd.xvg" "${PREFIX}_chain${c}_rmsd.png" \
        || echo "[WARN] chain ${c} RMSD failed"
      echo "[CMD] bash ${SCRIPT_DIR}/multichain_chain_rmsf.sh $REF $FIT $NDX $g ${PREFIX}_chain${c}_rmsf.xvg"
      bash "${SCRIPT_DIR}/multichain_chain_rmsf.sh" "$REF" "$FIT" "$NDX" "$g" "${PREFIX}_chain${c}_rmsf.xvg" \
        && gromd-plot-xvg "${PREFIX}_chain${c}_rmsf.xvg" "${PREFIX}_chain${c}_rmsf.png" \
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
          && gromd-plot-xvg "${out}.xvg" "${out}.png" \
          || echo "[WARN] inter-chain ${a}-${b} distance failed"
      done
    done
  else
    echo "[WARN] chain index produced <2 chains; skipping per-chain extras"
  fi
fi

echo "[OK] Post-analysis done → ${ANALYSIS_DIR}/"
