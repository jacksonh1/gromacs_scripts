# Analysis Scripts

GROMACS T-REMD post-processing utilities will live here.

## Scripts

### `process_trajectory.sh` — PBC correction
```
bash process_trajectory.sh OUTDIR [REP]
```
Fixes periodic boundary condition artifacts in the raw trajectory using `gmx trjconv -pbc mol -center -ur compact`. Output: `OUTDIR/analysis/remd_rep<REP>_pbc.xtc`.

**`-pbc nojump` is intentionally not used.** `nojump` compares consecutive frames to detect box-crossing events, but T-REMD coordinate exchanges cause discontinuous jumps between frames that `nojump` misinterprets. `-pbc mol` is a per-frame operation and is safe for REMD trajectories.

---

### `strip_and_align_trajectory.sh` — strip waters/ions + backbone alignment
```
bash strip_and_align_trajectory.sh OUTDIR [REP] [REF_GRO]
```
Requires `process_trajectory.sh` first. Fits each frame to a reference structure (default: first frame of the PBC-corrected trajectory) using backbone atoms (N, CA, C, O), and outputs only protein atoms — dropping waters and ions. Output: `OUTDIR/analysis/remd_rep<REP>_fit.xtc`.

The two operations are combined in one `gmx trjconv` call: `-fit rot+trans` with Backbone as the fitting group and Protein as the output group.

---

### `remd_acceptance.py` — exchange acceptance rates
```
python remd_acceptance.py OUTDIR [--rep REP] [--plot]
```
Parses the "Replica exchange statistics" block from the GROMACS log (pre-computed by GROMACS at the end of each run) and reports per-pair empirical acceptance rates, mean Metropolis probabilities, and exchange counts. Outputs a console table and `OUTDIR/analysis/remd_acceptance.csv`. Use `--plot` to also write `remd_acceptance.png`.

Target acceptance rates for T-REMD are 20–30% per adjacent pair. Rates above 50% indicate replicas are too closely spaced in temperature.

Each QC metric is a separate script so they can be modified independently.

---

## Typical workflow

```bash
bash analysis_scripts/process_trajectory.sh         OUTDIR    # 1. fix PBC
bash analysis_scripts/strip_and_align_trajectory.sh OUTDIR    # 2. strip + align
python analysis_scripts/remd_acceptance.py          OUTDIR    # 3. QC: acceptance rates
```

---

Planned tools:
- RMSD / Rg extractor for the 300 K replica trajectory
- Empirical transition matrix analysis (mixing quality between temperature slots)
- Round-trip counter (requires per-frame walker trajectory reconstruction from `Repl ex` lines)

## Key difference from Amber T-REMD

In GROMACS T-REMD, each replica runs at a **fixed temperature** and **coordinates** are exchanged between neighboring replicas. This means:

- `prod/rep000/remd.xtc` **is** the 300 K constant-temperature trajectory — use it directly for analysis.
- `prod/rep001/remd.xtc`, `rep002/remd.xtc`, etc. are the higher-temperature trajectories.
- **No demux step is needed** to obtain the constant-temperature ensemble.

`gmx demux` generates a trajectory that follows a specific *configuration* as it walks through temperature space — this is useful for visualizing the random walk of a molecule across replicas, but is not the correct tool for obtaining the thermodynamic ensemble at a given temperature.
