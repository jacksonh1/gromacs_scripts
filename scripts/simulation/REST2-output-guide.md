# GROMACS REST2 Output Guide

Replica exchange with **solute tempering** (REST2), via PLUMED Hamiltonian replica
exchange. Starts from a designed/folded input pose — never an unfolded state.
Primarily used for **stability/rigidity** (objective #1), **flexible-region**
identification (#2), and **variant comparison** (#3) — the same objectives as T-REMD.

**How it differs from T-REMD:** every replica runs at the *same physical temperature*
(`T_MIN`). Instead of heating the whole box, the **protein** force field is scaled by
`lambda_i = T_MIN / T_i` across a geometric **effective**-temperature ladder
(`T_MIN`…`T_MAX`). Water and ions are unscaled. Because only the solute is "heated,"
REST2 needs far fewer replicas than T-REMD for the same effective range.

**`rep000` (lambda = 1) is the true, unscaled `T_MIN` ensemble** — analyze it directly,
exactly like `rep000` in T-REMD. No demux needed.

Generated from `REST2-gromacs.sbatch`, which uses the **GROMACS 2023.5 + PLUMED** build
(`REST2_GMXRC` in `site_config.sh`) — the 2024.3 build's hrex is broken. The top-level
output directory is set by `OUTDIR` in the submit script, e.g.:

```
outputs/output_REST2/helix_fusion-20ns-REST2-300-450Keff-24reps-NVT-exf-1ps/
```

---

## Directory Tree

```
OUTDIR/
├── build/                        # Step 3  — System building                     [REAL]
├── em/                           # Step 4  — Energy minimization                  [REAL]
├── topol/                        # Step 6  — Per-replica SCALED topologies        [REAL]
│   ├── processed.top             #           grompp -pp expansion (unscaled)
│   ├── marked.top                #           protein atom types tagged "_" (mark_hot_region.py)
│   └── topol_rep000.top … NNN    #           plumed partial_tempering lambda_i  (rep000 = lambda 1)
├── analysis/                     # Post-analysis (rep000)                         [REAL]
├── logs/                         # mdrun stdout logs (incl. mdrun_preflight.log)  [REAL]
├── density/  ─▶ scratch          # Step 5  — NPT density equilibration (iterative)  [SYMLINK]
├── equil/    ─▶ scratch          # Steps 7–8 — Per-replica equilibration (all at T_MIN, scaled H)
│   ├── rep000/ … repNNN/
├── prod/     ─▶ scratch          # Steps 9–11 — REST2 production                  [SYMLINK]
│   ├── rep000/                   #           lambda=1 → the physical T_MIN ensemble
│   │   ├── rest2.tpr  rest2.log  rest2.gro  rest2.cpt  rest2.xtc
│   │   └── plumed.dat            #           minimal (hrex needs -plumed; no CVs)
│   └── rep001/ … repNNN/
├── parameters.txt                # Summary of all job parameters (incl. the lambda ladder)
└── <OUTBASE>_final_rep000.pdb    # Final structure from replica 000
```

**Output model (folder-symlink).** Small dirs (`build/ em/ topol/ analysis/ logs/`,
`parameters.txt`, final PDB) are **real** in `OUTDIR`; the bulk stage dirs (`density/
equil/ prod/`) are **folder symlinks into scratch** (`SCRATCH_DIR`), so mdrun writes
straight onto scratch and analysis reads `prod/rep000/rest2.xtc` through the symlink.
On success the run copies the real dirs (incl. `topol/`) into `SCRATCH_DIR` so the
scratch archive stands alone. Set **`SYMLINK_BULK=0`** to keep the stage dirs **real in
`OUTDIR`** with no scratch offload; everything else here is identical.

---

## Key files

| Path | What it is |
|------|------------|
| `prod/rep000/rest2.xtc` | **The analysis target** — the unscaled `T_MIN` ensemble (on scratch via the `prod/` symlink) |
| `prod/rep000/rest2.log` | Exchange log; the `Replica exchange statistics` block has per-pair acceptance |
| `topol/topol_rep*.top` | The scaled Hamiltonians (reproducible from `marked.top` + lambda) |
| `analysis/remd_acceptance.csv` | Per-pair exchange acceptance (target ~20–40%; retune the ladder if far off) |
| `parameters.txt` | Includes the lambda ladder and effective-T range actually used |

---

## Notes / gotchas

- **The engine fail-loud checks that hrex actually exchanges** (a ~20-attempt pre-flight,
  Step 10) *before* the real production, so a wrong/broken build cannot silently produce a
  fake REST2 run. If it aborts there, check `REST2_GMXRC` points at the 2023.5 build.
- **Acceptance too high (e.g. >60%)** → the ladder is over-dense; widen `T_MAX` or reduce
  `REPLICAS`. **Too low (<15%)** → add replicas or lower `T_MAX`. Tunable per system, like
  T-REMD's temperature ladder.
- **Runtime build**: REST2 loads `deprecated-modules` + `cuda/12.4.0` and prepends
  `HCOLL_COMPAT_DIR` to `LD_LIBRARY_PATH` (so `gmx_mpi` and the `plumed` CLI run on GPU
  nodes lacking `/opt/mellanox/hcoll`). All configured in `site_config.sh`.
- Production basenames are `rest2.*`; `run_analysis.sh` auto-detects REST2 from
  `prod/rep000/rest2.tpr` and treats `rep000` as the ensemble of interest.
