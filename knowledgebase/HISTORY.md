# History

Dated log of work. Newest entries on top. One entry per session/day.

---

## 2026-07-29 (latest) — REST2 engine written, tested end-to-end, DONE

- Wrote the full REST2 pipeline: `REST2-gromacs.sbatch` (15-step engine on the 2023.5 build,
  inline scaled-topology generation, constant-T equil/prod under scaled Hamiltonians, `-hrex`
  production, fail-loud acceptance preflight), `mark_hot_region.py`, `check_hrex_acceptance.py`,
  `submit_REST2.sh`, `REST2-output-guide.md`; extended `validate_params.py` (rest2),
  `run_analysis.sh` (REST2 mode), `remd_acceptance.py` (auto-detect rest2.log), `site_config.sh`
  (REST2_GMXRC/REST2_CUDA_MODULE/HCOLL_COMPAT_DIR). README REST2 → implemented.
- Defaults: 24 reps, 300→450 K effective (per user). hcoll compat libs staged at ~/opt/hcoll_compat.
- **End-to-end test PASSED** (job 19200028): all 15 steps, hrex preflight "7/7 pairs exchanging",
  production acceptance .27–.58, analysis wrote remd_acceptance.csv + RMSD/Rg/RMSF/DSSP.
- Bugs caught only by running (not `bash -n`): hrex-flag grep pattern, `TEMPS_LIST` `set -u`,
  hrex `pipefail`/SIGPIPE false-negative, cross-FS `OUTDIR` rename → all fixed; last two also
  new CLAUDE.md gotchas.

## 2026-07-29 (later) — REST2 UNBLOCKED: built GROMACS 2023.5+PLUMED, hrex verified working

- Built a dedicated **GROMACS 2023.5 + PLUMED 2.9.4** at `$HOME/opt/gromacs/2023.5-plumed`
  (2024.3 untouched for T-REMD/MD). Recipe: `scripts/installation/install_gromacs-2023.5-plumed.sh`.
- Cleared a build gauntlet (2023-era code on a 2026 cluster), each fix in the script +
  a CLAUDE.md gotcha: plumed CLI dead on compute nodes (patch on login), CMake 4.x rejects
  old googletest (load cmake/3.24.3), cuda 12.9/13.x too new (cuda/12.4.0 via deprecated-modules),
  build node (Sapphire Rapids) newer than Skylake run nodes (force -march=skylake-avx512).
- **hrex VERIFIED (job 19186450):** sanity 2×identical-λ → exchange P=1.0 (was 0.00 on 2024.3);
  8-rep fine ladder 300→350 K → 65-82% acceptance. Also re-ran the scale=0.5 all-atoms control
  on the 2023.5 binary — reproduces 2024.3 to the last digit.
- Net: all science + toolchain verified. **Next: write the REST2 engine** (REST2-gromacs.sbatch,
  mark_hot_region.py, validate_params rest2, submit_REST2.sh, run_analysis 3rd mode, output guide),
  pointing at the 2023.5 build + a fail-loud nonzero-acceptance self-check.

## 2026-07-29 — REST2 BLOCKED: `-hrex` silently broken on the 2024.3-plumed build

- Checked **acceptance** (not just completion) of the hrex smoke test: **0 exchanges,
  `dE_term ≡ 0`** even on a fine 8-rep/2%-λ/100 ps ladder. Retracts yesterday's
  premature "GPU-resident hrex is fine" (that only checked the job finished).
- Diagnosed: NOT a GPU/offload issue (0 exchanges on `-update cpu` and full `-nb cpu`
  too). Topologies correct, PLUMED engages. Root cause = the **PLUMED 2.9.x hrex patch
  for GROMACS 2024.x is an incomplete port** (`md.cpp` `GREX cacheLocalUSwap` energy
  re-eval wrong on the 2024 force API). Confirmed by two GROMACS-forum threads.
- **Blocker:** can't build REST2 on `2024.3-plumed_2.9.4`. Working combos: GROMACS
  2023.5 + PLUMED 2.9.x (safe) or GROMACS 2025 + PLUMED 2.10 native. Awaiting user's
  build decision.
- Added a **CLAUDE.md gotcha** (hrex silently broken here) + the earlier hcoll/rocky8
  gotcha. Plan doc + memory updated. T-REMD/MD unaffected.

## 2026-07-28 — REST2 pipeline: design + mechanism validation; KB scaffolded

- **REST2 pipeline design.** Broad plan agreed: a separate `REST2-gromacs.sbatch`
  engine, ~85% the T-REMD engine. Decisions locked (single-chain whole-protein hot
  region; `rest2.*` basename + a third `run_analysis.sh` mode). Full plan →
  `plans/REST2-pipeline.md`.
- **Confirmed the environment** does REST2: GROMACS 2024.3-plumed_2.9.4 exposes
  `-plumed`/`-hrex`/`-replex`; `plumed partial_tempering` present.
- **Validated the topology-scaling mechanism (scale=1.0 test PASSED):** drafted
  `mark_hot_region.py` (append `_` to the type col of the whole `Protein*`
  moleculetype), ran `partial_tempering 1.0` on the example `helix_fusion` build,
  and confirmed `mdrun -rerun` energies are identical to the bit vs the unscaled
  FF on amber99sb-ildn.
- **GPU `-hrex` runtime test submitted** (does it need `-update cpu`?) — running as
  job 19122812. First attempt died as a no-op because it was staged on login-local
  `/tmp` (compute nodes don't share it); restaged on shared `/orcd/data/keating`.
- **Knowledgebase scaffolded:** renamed the copied "FragForge" README, added
  `DECISIONS.md`, `SCIENCE.md`, and `plans/`. Gotchas stay in `CLAUDE.md`.
- **Open:** GPU hrex update-mode result; then write the engine + validator +
  submit script + output guide.
