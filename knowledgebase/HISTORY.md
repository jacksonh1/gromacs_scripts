# History

Dated log of work. Newest entries on top. One entry per session/day.

---

## 2026-09-02 (latest) — Barostat settings threaded through all three engines

Follow-on from the NPT-REMD discussion. Two conclusions from that discussion are recorded
here because they reverse earlier claims:

- **"NVT is the T-REMD convention" was wrong and is withdrawn.** Sugita & Okamoto 1999 was
  NVT, but practice is genuinely mixed and NPT is common. NVT also has its own defect: the
  box is frozen at the `T_MIN` equilibrium volume, so the top of the ladder runs over-dense
  at several hundred bar. The real risk in `ENSEMBLE=NPT` is not the ensemble, it is that
  TIP3P boils near 359-367 K at 1 bar while `T_MAX` defaults to 400 — the hot replicas sit
  ~40 K superheated. Deferred, not fixed; `REF_P` is now the lever if we want it.

- **`REF_P`/`TAU_P` were REMD/REST2-only params and are now shared** (`_SHARED_JOB`, checked
  for all three engines). They were also mislabelled "NPT only" everywhere. They are not:
  `density/` is pressure-coupled in every engine regardless of `ENSEMBLE`, since it is the
  stage that *sets* the box. `REF_P` is therefore the pressure the box production inherits
  was equilibrated at, whatever `ENSEMBLE` says. `ENSEMBLE` controls only the stages after
  `density/`: NPT keeps the barostat on in `equil/`+`prod/`, NVT turns it off and freezes the
  box at whatever `density/` ended on.

- **The density stage was ignoring both params.** All three engines hardcoded `tau-p = 1.0`
  / `ref-p = 1.0` in the density mdp (MD also in `relax/` and `prod/`), so setting `TAU_P=2`
  silently did nothing there and raising `REF_P` would have had `density/` equilibrate at
  1 bar and `equil/` re-do the work at the new pressure. Now `${TAU_P}`/`${REF_P}` throughout.

- **`TAU_P` default 1.0 -> 2.0.** C-rescale is thermodynamically consistent — unlike Berendsen
  it reproduces the correct NPT distribution at *any* `tau-p` — so this changes the volume's
  relaxation time, not the ensemble. The reason not to go below ~1 ps is that the pressure
  autocorrelation time is 0.1-0.5 ps, so a short `tau-p` makes the barostat chase instantaneous
  virial noise, and in `density/` (with `-DPOSRES` and `refcoord-scaling = com`) that does
  spurious work against the restraints. This also changes **MD production** `tau-p` 1.0 -> 2.0.

- MD gained `REF_P`/`TAU_P`, which it had never had. Done deliberately rather than leaving MD
  hardcoded: the planned validation is to measure the isothermal compressibility
  kappa_T = <dV^2>/(kB T <V>) on plain-MD-NPT at 300 K and on REMD-NPT `rep000` and compare,
  which is confounded if the two engines run different `tau-p`. Note the target is TIP3P's own
  value (~5.7-6.1e-5 /bar), not experiment (4.5e-5) — and the `compressibility` mdp line is
  only the barostat coupling scale, it does not constrain the measured value.

- Doc sweep alongside: the REMD stage list in `README.md` and `DECISIONS.md` still read
  `build/ em/ density/ equil/ prod/` (no `heat/`), and `README.md`/`GOTCHAS.md` still called
  `heat/` MD-only — all stale since the heat stage landed earlier the same day.

- Not done, by request: the `equil/` -> `density2/` rework. Filed in `knowledgebase/TODO.md`.

## 2026-09-02 — Equilibration audit: heat stage for REMD/REST2, density-loop fixes

Audit of the density + equilibration stages across all three engines. Five real defects, all
fixed; the NPT-REMD *physics* concerns (superheated water at `REF_P=1 bar` for a `T_MAX=400`
ladder, unverified per-replica box relaxation) are **deferred, not addressed** — pending a
conventions discussion.

- **`heat/` added to REMD and REST2** (new STEP 5), matching MD: NVT at `T_MIN`, restrained,
  velocities generated there. Previously both engines went `em/ → density/` directly, starting
  C-rescale on a minimized configuration whose instantaneous virial pressure is far from
  equilibrium; the `dt=0.001` first density segment was the crash guard for that. With heat in
  front, every density segment is now identical (uniform length, `DT_PS`), which also removed
  the dead `GEN_VEL`/`GEN_TEMP`/`SEG_STEPS`/`SEG_DT` branch — leftovers from the `SEED`
  retrofit, when `gen-vel`/`gen-temp` were folded into `CONT_FLAGS`.
- **Density segments now chain by checkpoint** (`grompp -t`) in REMD/REST2; only MD did. A
  `.gro` carries no thermostat/barostat state, so the barostat was restarting cold at every
  boundary and kicking the volume — the quantity the convergence test measures. New gotcha.
- **POSRES reference pinned to `em/em.gro`** in all three density loops (MD included). They
  passed `-r "$LAST_GRO"`, re-anchoring the restraint to the already-drifted structure each
  segment, so cumulative drift off the designed pose was never resisted. New gotcha.
- **Convergence is now a plateau/slope test** over the trailing `DENSITY_MIN_SEG` segments
  (`scripts/simulation/density_converged.py`, shared by all three engines, + `tests/`). The old
  consecutive-segment difference passes a sustained sub-tolerance drift on every comparison
  while the box contracts several percent. New gotcha.
- **`HEAT_NS` is now a shared param** (`_SHARED_JOB`), and `> 0` rather than `>= 0` — every
  engine's density stage resumes from `heat/heat.cpt`. `parameters.txt` records it, plus whether
  density ended at a plateau or the `DENSITY_MAX_SEG` cap.
- Renumbered REMD STEPs 5–14 → 6–15 and REST2 5–16 → 6–17 (incl. `PRESERVE_FROM_STEP`), updated
  both output guides, `PARAMETERS.md`, `SCIENCE.md`, `CLAUDE.md`, the config examples, and
  corrected a stale claim in `../archive/amber_vs_gromacs_pipeline_comparison.md` (REMD restraints are
  released at the *start* of `equil/`, not at production start — the equil mdp has no `-DPOSRES`).
- **Validated on the cluster.** Short smoke jobs on `helix_fusion` (MD 21833585, REMD 21833586)
  both `COMPLETED 0:0` through every step incl. analysis and the scratch finalize. Confirmed:
  `heat/` is created as a scratch symlink and writes `heat.gro`/`heat.cpt`; density `grompp` logs
  show "Reading Coordinates, Velocities and Box size from old trajectory" (the `-t` checkpoint);
  `parameters.txt` shows `Heat ns` and `Density convergence segments: 8 (plateau reached)`.
- **The POSRES fix is visible in the logs.** `grompp`'s "center of mass of the position restraint
  coord's" is now *constant* across all 8 density segments (3.808 4.124 1.940), where the archived
  old run walked it every segment (4.113 → 4.075 → 4.057 → 4.081 → 4.072 in x). Caveat: part of
  that old movement is isotropic box rescaling under `refcoord-scaling = com`, so this confirms the
  code change took effect but does not by itself measure protein drift.
- **The plateau test agrees with the old test on a well-behaved box** — both stop at segment 8 —
  so it is not spuriously stricter. MD drift 0.16 %, REMD 0.10 %, tol 0.50 %.
- Follow-up caught by reading the job logs: the drift diagnostic was going to stderr while the
  per-segment volumes went to stdout, splitting one decision across the SLURM `.out` and `.err`
  files. `density_converged.py` now prints `<verdict> <reason>` on a single stdout line and each
  engine splits it, so the reason lands next to the volumes. 53 tests pass.
- A rerun is NOT expected to reproduce an old run bit-for-bit (random seeds, GPU reduction order);
  the check here is that the new stages run and produce the right shape of output, not a diff.
- **Coverage pushed past the one small system**, after the first round only exercised `helix_fusion`:
  the `DENSITY_MAX_SEG` cap branch (forced, 21834958 — warning fires, job still `COMPLETED 0:0`,
  `parameters.txt` annotated); `1a22-fixed` at 949 021 atoms / 2 chains / 9580 nm^3 (21834959);
  CHARMM36m through the new REMD `heat.mdp` (21834960 — force-switch block correct). Heat handled
  the large strained multi-chain system with **zero** LINCS warnings. The plateau test converged at
  segment 8 across a 90x range of box volume (115 -> 9580 nm^3). A Monte Carlo over box sizes puts
  the noise-only false-negative rate at 0.03 % worst case (smallest box), falling toward zero as
  boxes grow — but the test is ~7x more sensitive than the old one, so expect some systems to run
  past 8 segments now.
- Ladder note from the big run: `1a22-fixed` gave 0 % exchange acceptance on an 8-replica
  300-350 K ladder — correct physics, not a defect (replica count must scale as sqrt(N)); the
  archived 48-replica production run shows a healthy uniform .33-.39 across all 47 pairs. Read
  acceptance out of the log's `Replica exchange statistics` block; do not re-derive it.

## 2026-09-02 — Two latent defects found while smoke-testing the above

- **`step*.pdb` crash dumps were landing in the SUBMIT DIRECTORY.** GROMACS writes constraint-failure
  dumps to the cwd, which for an sbatch is where you ran `sbatch` — not `OUTDIR`. On `1a22-fixed`
  (949 021 atoms) each is 75 MB; one run dropped 131 MB next to the submit scripts, on the very
  tight-quota pool the `SYMLINK_BULK` model exists to protect. Fixed: each engine now does
  `cd "${SCRATCH_DIR:-$OUTDIR}"` **immediately before EM**. Verified (21836304): 0 files in the
  submit dir, all 5 dump pairs (750 MB) on scratch. **Placement is load-bearing** — putting the `cd`
  at the obvious spot (right after the stage dirs are made) moves cwd off the pool while
  `solvate`/`genion` are still running and their cross-filesystem topology `rename()` kills the
  build 8 s in (measured, 21836215). New gotcha; it links to the existing OUTDIR-filesystem one.
  Not caused by the heat-stage work: the dumps come from **EM** on strained input, and the archived
  MD run on the same structure did it too.
- **`FORCE` was a documented parameter that had never done anything — removed.** Declared in the
  Amber REMD scripts (`amber_REMD/amber_scripts/REMD-N*T-singlenode.sbatch`), copy-inherited into
  the 2026-03-25 GROMACS port, then into this engine, then into REST2 (written by copying REMD).
  All 15 copies across the workspace contain exactly one occurrence: the self-referential
  `FORCE="${FORCE:-0}"` default. Nothing ever read it. The Amber REAF/REST2 scripts have it
  *commented out*, so it was spotted once and never cleaned from the REMD line. `validate_params.py`
  later gave it a real `flag01` check and `PARAMETERS.md` described it as "Overwrite an existing
  OUTDIR" — a behavior invented to fit a dead variable, which is why `FORCE=1` does not rescue a
  re-run into an existing `OUTDIR`. Removed from both engines, the validator allowlist, and the docs.
  **Still open:** a failed run leaves stale stage symlinks in `OUTDIR` and the next run dies at
  `ln -s ... File exists`. There is now no override for that; clear the `OUTDIR` by hand.
- **Made the dumps discoverable without weakening `PRESERVE_SCRATCH_FROM`.** Moving them to
  scratch fixed the quota problem but created a visibility one (nobody browses scratch), so
  `report_crash_dumps.sh` now writes a few-KB `OUTDIR/CONSTRAINT_FAILURES.txt` — count, which
  stage's log recorded them, full paths, how to diff a b/c pair — plus a banner in the job log.
  It clears a stale report when a run is clean, and always exits 0.
  **A first cut had the ERR trap force-preserve scratch whenever dumps existed; that was wrong
  and was reverted.** `never` has to mean never, and since a strained structure dumps during EM
  on *every* run (`1a22-fixed` does), the "override" would have fired constantly — disabling the
  knob for that whole class of input and stranding the entire multi-GB `density/equil/prod` tree,
  not just the dumps. Final split: the policy governs the bulky PDBs; the cheap report always
  survives, and when scratch is cleaned the report says so and names the setting
  (`PRESERVE_SCRATCH_FROM=always`) that would have kept them. Verified end-to-end on
  `1a22-fixed` (21837272, `COMPLETED 0:0`): banner in the log, breadcrumb in OUTDIR, all 5 events
  correctly attributed to `em/em.log`, 326 MB of PDBs left on scratch.

---

## 2026-07-29 — REST2 engine written, tested end-to-end, DONE

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
