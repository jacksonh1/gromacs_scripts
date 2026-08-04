# Plan: REST2 pipeline

**Status:** DONE (2026-07-29). Engine written and verified end-to-end. REST2 is the third
engine alongside T-REMD/MD, on a dedicated GROMACS 2023.5 + PLUMED 2.9.4 build (hrex works
there; the 2024.3 build's hrex is silently broken). Plain T-REMD/MD stay on 2024.3-plumed.

**Delivered files:** `scripts/simulation/REST2-gromacs.sbatch` (engine, 15 steps),
`mark_hot_region.py` (protein marker), `check_hrex_acceptance.py` (fail-loud gate),
`REST2-output-guide.md`; `example/submit_jobs/submit_REST2.sh`; `validate_params.py`
(`--engine rest2`); `run_analysis.sh` (REST2 mode) + `remd_acceptance.py` (auto-detects
rest2.log); `site_config.sh` (REST2_GMXRC / REST2_CUDA_MODULE / HCOLL_COMPAT_DIR);
`install_gromacs-2023.5-plumed.sh`; README flipped REST2→implemented.

**End-to-end test PASSED** (job 19200028, 8 reps × 50 ps, mit_normal_gpu): all 15 steps,
inline `partial_tempering` on the compute node (hcoll LD_LIBRARY_PATH), hrex preflight gate
fired ("7/7 pairs exchanging"), production acceptance .27–.58 over 300→450 K, analysis
auto-detected REST2 and wrote remd_acceptance.csv + RMSD/Rg/RMSF/DSSP. (Clustering warned —
6-frame test too short; non-fatal, fine on real runs.)

REST2 = Replica Exchange with Solute Tempering (v2), via the PLUMED-patched
GROMACS. It is ~85% the existing T-REMD engine; the science is in `../SCIENCE.md`,
the design rationale in `../DECISIONS.md`.

## Locked decisions

- Used exactly like REMD: **single-chain, whole-protein hot region.** No complexes,
  no residue-subset selection.
- Production basename **`rest2.*`** + a **third detected mode** in `run_analysis.sh`
  (rep000/λ=1 analysis is identical to REMD's).
- A **separate parallel engine** `REST2-gromacs.sbatch`, not a branch in REMD.

## Environment (confirmed present)

- GROMACS **2024.3-plumed_2.9.4** at `~/opt/gromacs/2024.3-plumed` — `gmx_mpi mdrun`
  exposes `-plumed`, `-hrex`, `-replex`.
- PLUMED **2.9.4** on PATH; `~/plumed.sh` sets `PLUMED_KERNEL`/`LD_LIBRARY_PATH`.
- `site_config.sh` already carries `GMXRC` (plumed build) + `PLUMED_SH`.

## The mechanism (confirmed by the partial_tempering source + a scale=1.0 test)

`plumed partial_tempering <λ> < processed.top` reads a `grompp -pp` topology in
which "hot" atoms are flagged by **appending `_` to the atom-type column (col 2)**
of a `[ atoms ]` section, and scales charge ×√λ, LJ ε ×λ, 1-4 pairs, and proper
dihedrals accordingly. So the new step is:

```
grompp -pp processed.top          # once, unscaled, NO -DPOSRES (scaled tops carry no restraints)
mark_hot_region.py processed.top  # append "_" to type col of every Protein* [atoms] line
for i in replicas:
    plumed partial_tempering λ_i < marked.top > topol_i.top   # λ_i = T_MIN/T_i, λ_0 = 1.0
```

Then per-replica `grompp -p topol_i.top`, all replicas at `ref-t = T_MIN`, and
production `mdrun ... -plumed plumed.dat -hrex -replex`.

Constraints baked in:
- `-replex` steps must be a multiple of `nstlist` (hrex swaps only at neighbor-search
  steps; REMD's 500 % 10 already satisfies this).
- `REPLEX_PS >= 1 ps` gotcha applies (likely sharper under hrex).

## New engine = REMD engine + these deltas

| Step | Change vs REMD-gromacs.sbatch |
|---|---|
| 0 Env | also `source PLUMED_SH`; verify `plumed` + that mdrun accepts `-hrex` (fail loud) |
| 1 Params | reinterpret `T_MAX` = effective max solute T; lower default `REPLICAS` (~8–16) |
| 2 Ladder | same geometric math + emit λ_i = T_MIN/T_i; all `ref-t` = T_MIN |
| 3–5 Build/EM/Density | identical (unscaled topology is correct for density equil) |
| **NEW Topol** | `grompp -pp` → `mark_hot_region.py` → per-replica `partial_tempering` |
| 6/7 Equil | per-replica grompp uses `topol_i`; `ref-t` = T_MIN |
| 8/9 Prod | per-replica grompp uses `topol_i`; mdrun gains `-plumed plumed.dat -hrex` |
| 10/11 | identical (basename `rest2.*`) |
| 12 Analysis | `run_analysis.sh` third mode on `prod/rep000/rest2.tpr` |

## Files to add / touch

- `scripts/simulation/REST2-gromacs.sbatch` — new engine (clone of REMD)
- `scripts/simulation/mark_hot_region.py` — new; stdlib, fail-loud (crash on 0 marked)
- `scripts/simulation/validate_params.py` — add `--engine rest2` + allowlist
- `scripts/simulation/config_example.sh` — REST2 example
- `example/submit_jobs/submit_REST2.sh` — new
- `scripts/simulation/REST2-output-guide.md` — new user-facing guide
- `scripts/analysis/run_analysis.sh` — detect `rest2.*` (third mode)
- `site_config.sh` — already ready (source PLUMED_SH in the engine)

## Validation log

- **2026-07-28 — scale=1.0 correctness: PASSED.** Reusing the example `helix_fusion`
  build (protein moleculetype `Protein_chain_A`, 433 atoms marked), `mdrun -rerun`
  energies (Potential / LJ-SR / Coulomb-SR / Coulomb-recip) were identical to the
  bit vs the unscaled force field on amber99sb-ildn. Confirms `mark_hot_region.py`
  and `partial_tempering` are correct on the default FF. Draft script + drivers in
  the session scratchpad.

- **2026-07-28 — GPU `-hrex` runtime test: RESOLVED (job 19123138, node3006 L40S).**
  `-hrex` runs cleanly with the **default GPU-resident update** — 4/4 replicas finished
  in 8 s, no hang; `-update cpu` also works. **So GPU-resident is fine; `-update cpu`
  is NOT required.** The `Replica exchange statistics` block is written and parseable.
  - **Caveat — sampling NOT yet validated (RESOLVED 2026-07-29, see below):** at the
    time of this entry acceptance was 0 exchanges / ~0.00 prob across all 3 pairs (19
    attempts) — but that was on the silently-broken 2024.3 build. On the fixed GROMACS
    2023.5 + PLUMED 2.9.4 build acceptance is 0.65–0.82 on a fine ladder (and P=1.0 on
    the identical-λ sanity control), so sampling is now validated. Historical note kept
    for the chronology; it is not the current status.
  - **Two environment findings for the real engine (not just the smoke test):**
    1. **hcoll:** `gmx_mpi` has a baked-in `NEEDED` on `libhcoll.so.1`/`libocoms.so.0`
       (from `/opt/mellanox/hcoll/lib`, inert — MPI routes via `libmpi`). Some
       `mit_normal_gpu` nodes (e.g. node4207) lack it → `error while loading shared
       libraries`. Fix: stage both libs on shared storage + prepend `LD_LIBRARY_PATH`.
       Deps are all standard RDMA/system libs. **→ CLAUDE.md gotcha candidate.**
    2. **OS:** binary is Rocky-8-built; add `-C rocky8` (all 67 partition nodes are
       rocky8 today, so belt-and-suspenders/future-proofing).
  - Aborted earlier attempts: 19122222 (login-`/tmp` not shared w/ compute → no-op in
    4 s; restaged on shared Keating scratch); 19122812 (same, before restage).

- **2026-07-29 — scale=0.5 all-atoms scaling test PASSED.** Marked ALL atoms
  (test-only `mark_all_atoms.py`), `partial_tempering 0.5`, `mdrun -rerun` energies:
  Bond & Angle ratio 1.000 (unchanged); Proper-Dih., Per.-Imp.-Dih., LJ-14,
  Coulomb-14, LJ-SR, Disper.-corr., Coulomb-SR, Coul.-recip. all ratio 0.500;
  Potential 0.4965 (sum, bonds/angles stay full). Matches the partial_tempering
  help exactly → topology scaling is correct in BOTH directions (×1.0 identity,
  ×0.5 halving). The scaling was never the issue; only the 2024 hrex patch is.

- **2026-07-29 — hrex VERIFIED WORKING on the new GROMACS 2023.5 + PLUMED 2.9.4 build
  (job 19186450).** Built to `$HOME/opt/gromacs/2023.5-plumed` (cuda/12.4.0,
  cmake/3.24.3, forced `-march=skylake-avx512`, source pre-patched on login node).
  Sanity (2 identical-λ replicas, CPU): average exchange probability **1.0** — vs 0.00
  on 2024.3 for the same test. Fine ladder (8 reps, 300→350 K eff, GPU): per-pair
  acceptance **0.65-0.82** (over-dense — production should space replicas out toward
  ~30-40%). scale=0.5 all-atoms control also reproduced on 2023.5 to the last digit.
  Build gauntlet (all fixed in `install_gromacs-2023.5-plumed.sh`): plumed CLI dead on
  compute nodes → patch on login; CMake 4.x rejects old googletest → cmake/3.24.3;
  cuda 12.9/13.x too new → cuda/12.4.0 via deprecated-modules; build CPU (Sapphire
  Rapids) newer than run nodes → forced Skylake-AVX512.

- **2026-07-29 — `-hrex` is SILENTLY BROKEN on the 2024.3 build. (root cause, now resolved)** After the
  smoke test "ran", checking *acceptance* (not just completion) showed **0 exchanges /
  `dE_term ≡ 0.000e+00`** for every pair, even on a fine ladder (8 reps, 2% λ steps,
  100 ps, job 19123560). Diagnosis:
  - NOT a GPU issue — `-update cpu` AND full `-nb cpu` both give 0 exchanges (job 19172350).
  - Topologies are correct (charges scale as √λ) and PLUMED engages ("GROMACS-like
    replica exchange is on"). GROMACS logs "Replica exchange in temperature" at 8×300 K
    and, per `replicaexchange.cpp`, zeroes its own delta (`if(plumed_hrex) delta=0.0`)
    expecting PLUMED to supply the Hamiltonian delta via `md.cpp`'s `GREX cacheLocalUSwap`
    energy re-evaluation. On the **PLUMED 2.9.4 → GROMACS 2024.3** patch that re-eval
    returns a wrong value (incomplete 2024 port; the patch is structurally near-identical
    to the working 2023.5 one, only cosmetic member renames), so no real delta is added.
  - **Corroborated + authoritative root cause:** PLUMED dev **Giovanni Bussi** on the
    plumed-users list ("REST2 simulations not exchanging"): **"hrex is not yet supported
    in the 2024 patch."** The 2024.3 patch supports general PLUMED use but NOT the hrex
    energy machinery. PLUMED GitHub issue #1326 (GROMACS 2024.2 + PLUMED 2.9.2, `-hrex`
    → `dplumed=2.6e+05` kT → P=0; without `-hrex` → dplumed=0 → P=1) is OPEN, no fix.
    Reporting note: GROMACS zeroes its own `dE_term` under hrex by design; the real
    broken quantity is PLUMED's `dplumed` (huge). We only grepped `dE_term` (=0) — same
    conclusion. **Only confirmed-working combo: GROMACS 2023.5 + PLUMED 2.9.x** (scale-1.0
    sanity → P=1.0 there). GROMACS 2025 + PLUMED 2.10 native interface: unverified for hrex
    (2.10 changelog advertises multiple-walkers, not hrex). Per Bussi, 2023 is the newest
    GROMACS with working hrex — not an arbitrary downgrade.
  - → **CLAUDE.md gotcha added.** REST2 engine must (a) be built against a working
    GROMACS+PLUMED, and (b) gate on a nonzero-acceptance self-check (fail loud).

## Resolution — all complete (engine shipped)

Every item below was resolved; REST2 is a production engine, not a work in progress. Kept
as a record of what was closed out.

0. **Build fixed (was the blocker).** Chose Option A — a dedicated **GROMACS 2023.5 +
   PLUMED 2.9.4** build at `$HOME/opt/gromacs/2023.5-plumed`, alongside the existing 2024.3
   (kept for T-REMD/MD). `install_gromacs-2023.5-plumed.sh` shipped; builds against
   **cuda/12.4.0** (`deprecated-modules`), and the engine loads the same via `REST2_CUDA_MODULE`.
   hrex verified working (job 19186450): identical-λ sanity P=1.0, fine-ladder acceptance
   0.65–0.82. (Option B, GROMACS 2025 + PLUMED 2.10, not pursued — would migrate the whole
   toolchain.)
1. ~~GPU `-hrex` update-mode~~ — resolved: GPU-resident update works (no `-update cpu` needed).
2. **Acceptance validated** on the fixed build (nonzero exchanges on the fine ladder).
3. **λ ladder + replica count** — geometric effective-T ladder implemented; per-system
   acceptance tuning is inherent usage (like REMD's temperature spacing), not a gap.
4. **Charge-scaling / PME neutral background** — standard REST2 practice; noted in
   `REST2-output-guide.md`.
5. **Engine carries** the hcoll `LD_LIBRARY_PATH` handling, `-C rocky8`, `-dlb no`, and the
   fail-loud `check_hrex_acceptance.py` nonzero-acceptance gate (a broken build can never
   again masquerade as a working REST2 run).
6. **Engine written** — `scripts/simulation/REST2-gromacs.sbatch` (15 steps).
