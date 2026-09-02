#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# report_crash_dumps.sh — surface GROMACS constraint-failure dumps where they
#                         will actually be seen
# ─────────────────────────────────────────────────────────────────────────────
# GROMACS writes step<N>b.pdb / step<N>c.pdb (coordinates before and after
# constraining) into the process cwd whenever LINCS/SETTLE reports an excessive
# deviation. The engines point that cwd at SCRATCH_DIR so the files — 75 MB each
# on a large system — stay off the tight-quota pool.
#
# The cost of that is discoverability: nobody browses the scratch volume, so a
# run that quietly threw constraint failures and then finished would look clean.
# This script closes that gap. The bulky PDBs stay on scratch; a small, loud
# breadcrumb goes into OUTDIR, which is the directory people actually open.
#
# Usage:
#   bash report_crash_dumps.sh DUMP_DIR OUTDIR [retained|deleted]
#
# The third argument says what is about to happen to the dumps themselves, so the
# report can be honest about whether the paths it lists will still exist. The ERR
# trap passes "deleted" when PRESERVE_SCRATCH_FROM says to clean scratch: the small
# report still lands in OUTDIR (a few KB), but the 75 MB PDBs go, per policy.
#
# Writes OUTDIR/CONSTRAINT_FAILURES.txt and prints a banner to stdout (the job
# log) when dumps exist. Silent when there are none, and it REMOVES a stale
# report from an earlier run so the file's presence always means "this run".
# Always exits 0: this is diagnostics, and it must never fail a good job.
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

DUMP_DIR="${1:?Usage: report_crash_dumps.sh DUMP_DIR OUTDIR [retained|deleted]}"
OUTDIR="${2:?Usage: report_crash_dumps.sh DUMP_DIR OUTDIR [retained|deleted]}"
FATE="${3:-retained}"
REPORT="${OUTDIR%/}/CONSTRAINT_FAILURES.txt"

shopt -s nullglob
dumps=( "${DUMP_DIR%/}"/step*.pdb )
shopt -u nullglob

if (( ${#dumps[@]} == 0 )); then
  rm -f "$REPORT"          # a previous run's report must not linger
  exit 0
fi

# Each failure writes a PAIR (b = before constraining, c = after), so events =
# files/2. Report events, since that is what the GROMACS log counts.
events=$(( ${#dumps[@]} / 2 ))
total_h=$(du -ch "${dumps[@]}" 2>/dev/null | tail -1 | cut -f1)

{
  echo "GROMACS CONSTRAINT FAILURES DETECTED"
  echo "===================================="
  echo
  echo "This run triggered ${events} constraint-failure dump(s) (${#dumps[@]} files, ${total_h})."
  echo "GROMACS writes these when LINCS/SETTLE reports an excessive deviation —"
  echo "i.e. something moved much further in one step than it should have."
  echo
  echo "WHAT IT USUALLY MEANS"
  echo "  Most often: strained geometry in the input structure, hit during energy"
  echo "  minimisation. EM normally recovers and the run continues correctly — a"
  echo "  few dumps in em/ are common and not automatically a problem."
  echo "  Take it seriously if the dumps come from heat/, density/ or production,"
  echo "  or if the count is large: that is a simulation going unstable."
  echo
  echo "WHICH STAGE"
  for log in "${OUTDIR%/}"/em/em.log "${OUTDIR%/}"/heat/heat.log \
             "${OUTDIR%/}"/density/density_seg*.log "${OUTDIR%/}"/prod/*.log \
             "${OUTDIR%/}"/prod/rep000/*.log; do
    [[ -f "$log" ]] || continue
    n=$(grep -c "Wrote pdb files" "$log" 2>/dev/null || true)
    (( n > 0 )) && printf "  %-4s in %s\n" "$n" "$log"
  done
  echo
  if [[ "$FATE" == "deleted" ]]; then
    echo "THE FILES — ALREADY DELETED"
    echo "  The dumps were on scratch and have been removed with it, because"
    echo "  PRESERVE_SCRATCH_FROM said to clean scratch for a failure at this stage."
    echo "  This report is kept so the failure is still on record. To retain the PDBs"
    echo "  themselves on the next run, submit with PRESERVE_SCRATCH_FROM=always."
    echo
    echo "  They were:"
  else
    echo "THE FILES (on scratch — NOT backed up, and purged when the scratch dir is)"
  fi
  for f in "${dumps[@]}"; do
    printf "  %10s  %s\n" "$(du -h "$f" 2>/dev/null | cut -f1)" "$f"
  done
  echo
  echo "HOW TO INVESTIGATE"
  if [[ "$FATE" == "deleted" ]]; then
    echo "  The files are gone, so re-run with PRESERVE_SCRATCH_FROM=always to capture"
    echo "  them. Then diff a b/c pair (b = coordinates going in, c = after"
    echo "  constraining) to find the atoms that moved."
  else
    echo "  The b/c pair brackets one failure: b = coordinates going in, c = after"
    echo "  constraining. Diff a pair to find the atoms that moved:"
    echo "    diff ${dumps[0]} \\"
    echo "         ${dumps[1]:-<the matching c file>} | head"
    echo "  Then look at that residue in the input structure."
  fi
} > "$REPORT"

echo
echo "################################################################################"
echo "# [WARNING] ${events} GROMACS constraint-failure dump(s) — ${#dumps[@]} files, ${total_h}"
echo "#"
echo "# The run may still be fine (EM on strained input does this routinely), but"
echo "# something exceeded a LINCS/SETTLE tolerance. Details:"
echo "#   ${REPORT}"
echo "# Dumps are on scratch: ${DUMP_DIR%/}"
echo "################################################################################"
echo

exit 0
