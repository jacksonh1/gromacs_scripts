# GROMACS Structure-Characterization Pipelines - mostly stable

# warning - This pipeline is under active development


# upcoming features (not yet implemented):

- refactoring of the analysis layer to a python package
- implicit solvent support (currently only explicit water)
- adjustable force-field support (currently only AMBER99SB-ILDN)



GROMACS 2024.3 pipelines for characterizing **folded protein structures** on single-node GPU clusters (SLURM). The input can be any folded pose — a de novo design, a crystal/cryo-EM structure, a predicted model, or a mutant variant. Built for the Keating lab at MIT; configurable for any cluster via `site_config.sh`.

Three engines are provided:
- **T-REMD** (`submit_REMD.sh`) — temperature replica-exchange enhanced sampling
- **REST2** (`submit_REST2.sh`) — replica exchange with solute tempering (PLUMED Hamiltonian
  replica exchange): all replicas at one physical temperature, only the protein's force field
  scaled across an effective-temperature ladder. Fewer replicas than T-REMD for the same range.
- **Plain production MD** (`submit_MD.sh`) — single-temperature NPT production

## Purpose

The analysis treats the **input pose as the reference** — observables measure how far the structure moves from it, not whether it can be folded from scratch. The tools serve four objectives:

1. **Stability/rigidity** characterization of a given starting structure
2. Identifying **flexible regions** (per-residue)
3. **Variant comparison** — run the same protocol on several variants (mutants, designs, homologues, …) and see which best retains its starting conformation
4. **Bound-state sampling** — simulate a complex in its bound pose and sample the bound ensemble

Plain MD is mainly for **#4** (and optionally #1, #2); T-REMD is primarily for **#1–3**. The input pose is the reference the analysis is measured against (RMSD = drift from the input structure, RMSF = local flexibility).

---

## Prerequisites

- GROMACS 2024.3 compiled with MPI + CUDA, plus PLUMED — **full build instructions in [`scripts/installation/README.md`](scripts/installation/README.md)**
- A conda env for post-analysis (matplotlib, mdanalysis, numpy, …), created from `scripts/installation/environment.yml`, with the `gromd_analysis` package installed into it — see One-Time Setup. The sbatch engines use the system `python3` (standard library only) for inline temperature-ladder and convergence calculations, and activate the conda env only for the post-analysis/plotting steps.
- SLURM with GPU access

T-REMD and plain MD run on the GROMACS 2024.3 build. **REST2 needs a separate GROMACS 2023.5 + PLUMED build** — the PLUMED hrex patch is broken on GROMACS 2024 (runs but silently gives zero exchanges), so REST2 uses its own build (`REST2_GMXRC` in `site_config.sh`; see `scripts/installation/README.md` Step 5). Skip that build entirely if you don't run REST2. See `CLAUDE.md`.

---

## One-Time Setup

0. **Build GROMACS + PLUMED** (skip if your cluster already has a suitable build):
   see **[`scripts/installation/README.md`](scripts/installation/README.md)** — the ordered
   recipe for PLUMED, the 2024.3 build (T-REMD/MD), the 2023.5 build (REST2 only),
   the optional CHARMM36m force field, and how to port the build scripts to another
   cluster.

1. **Edit `site_config.sh`** in the repo root. Set:
   - `GMXRC` — path to your `bin/GMXRC` from the GROMACS installation
   - `SCRATCH_ROOT` — fast scratch storage root for trajectory files (~100 GB per job)
   - `CUDA_MODULE`, `OPENMPI_MODULE` — your cluster's module names
   - `CONDA_MODULE`, `GROMD_ENV` — module providing conda and the analysis env name

   (`REST2_GMXRC`, `GMXLIB`, `PLUMED_SH`, `HCOLL_COMPAT_DIR` matter only for REST2 /
   CHARMM36m — see the installation README.)

2. **Create the analysis conda env** (once per cluster, on a login node):
   ```bash
   bash scripts/installation/install_python_env.sh
   ```
   This builds `groMD_env` from `scripts/installation/environment.yml` **and**
   installs the `gromd_analysis` package into it (editable), which is what puts
   the `gromd-*` analysis commands on `PATH`. The pipeline activates the env
   automatically for the post-analysis steps.

   If the env already exists, just add the package:
   ```bash
   conda activate groMD_env && pip install --no-deps -e .
   ```

3. **Check the install** (optional, seconds):
   ```bash
   conda activate groMD_env && pytest
   ```

4. That's it. The pipeline scripts source `site_config.sh` automatically.

---

## Running a Job

1. Copy the example submit script for the engine you want and edit your parameters directly
   in it:
   ```bash
   # T-REMD:
   cp example/submit_jobs/submit_REMD.sh my_job.sh
   # edit my_job.sh: set PDB_IN, REPLICAS, T_MAX, TOTAL_NS, ENSEMBLE, FF, etc.

   # Plain MD:
   cp example/submit_jobs/submit_MD.sh my_job.sh
   # edit my_job.sh: set PDB_IN, T_SIM, TOTAL_NS, TRAJ_PS, RELAX_NS, FF, etc.

   # REST2:
   cp example/submit_jobs/submit_REST2.sh my_job.sh
   ```

2. Run it:
   ```bash
   bash my_job.sh
   ```

### Parameters

**→ [`docs/PARAMETERS.md`](docs/PARAMETERS.md) is the full reference** — every job
parameter for all three engines, with defaults, ranges, and which engine accepts it.

There are **two ways** to set job parameters, and they can be mixed:

- **In the submit script** (the common case) — edit the assignments at the top of your
  copy; the script passes them to `sbatch --export`.
- **In a separate config file** — all three engines take a config file as their **first
  argument** and source it. Start from `scripts/simulation/config_example.sh` and submit
  with `example/submit_jobs/submit_REMD_with_config.sh my_config.sh`, or call an engine
  directly: `sbatch -n 48 scripts/simulation/REMD-gromacs.sbatch my_config.sh`.
  A value in the config file **overrides** the same name passed via `--export`.

Anything set in neither place falls back to the engine default listed in
`docs/PARAMETERS.md`.

Cluster-level settings (GROMACS paths, modules, scratch root, force-field aliases) are
**not** job parameters — they live in `site_config.sh` and are edited once. SLURM
resources come from the `#SBATCH` headers in the engine script.

Parameters are validated at job start: a missing `PDB_IN`, an out-of-range or
non-numeric value, or (in a config file) a misspelled key fails the job immediately with
a clear `[ERROR]` message instead of silently falling back to a default. Typo detection
is **config-file only** — under `--export=ALL` a job inherits the whole environment, so
there is no manifest to check a misspelled name against.

Two parameter gotchas worth knowing before your first run:

- **`REPLEX_PS`** — never set below `1.0`. Sub-ps exchange deadlocks GPU-resident REMD
  mid-run with clean physics. The default is `1.0`; there is no reason to lower it.
- **`FF` under REST2** — CHARMM force fields are **rejected** by the REST2 engine; see
  [Force-field compatibility](#force-field-compatibility) below.

---

## Pipeline Overview

Both engines share the same system-building, equilibration philosophy (position restraints
held through equilibration to preserve the input pose), scratch handling, and post-analysis.
They differ in how production is sampled: T-REMD runs many temperature replicas with
exchanges; plain MD runs a single trajectory.

### T-REMD — `scripts/simulation/REMD-gromacs.sbatch`

| Step | Description |
|------|-------------|
| 0 | Load environment (modules, GROMACS) |
| 1 | Set parameters, create scratch and output directories |
| 2 | Compute geometric temperature ladder |
| 3 | Build system: pdb2gmx → editconf → solvate → genion |
| 4 | Energy minimization (steepest descent) |
| 5 | NPT density equilibration (iterative, convergence-checked) |
| 6 | Prepare per-replica equilibration inputs (NVT, or NPT if `ENSEMBLE=NPT`) |
| 7 | Run per-replica equilibration (all replicas in parallel via MPI) |
| 8 | Prepare REMD production inputs |
| 9 | Run T-REMD production (NVT or NPT) |
| 10 | Finalize outputs, create trajectory symlinks |
| 11 | Write parameters log |
| 12 | Post-analysis: acceptance rates + PBC/strip/align + RMSD/Rg/RMSF/DSSP + clustering (rep000) |

Stage folders: `build/ → em/ → heat/ → density/ → equil/ → prod/`. See
`scripts/simulation/REMD-output-guide.md` for a full description of all output files.

### Plain MD — `scripts/simulation/MD-gromacs.sbatch`

Single-temperature production, always NPT. Mainly for bound-state ensemble sampling (#4),
and optionally single-structure stability (#1) and flexible-region (#2) characterization.

| Step | Description |
|------|-------------|
| 0 | Load environment (modules, GROMACS) |
| 1 | Set parameters, create scratch and output directories |
| 2 | Build system: pdb2gmx → editconf → solvate → genion |
| 3 | Energy minimization (steepest descent) |
| 4 | Heat to `T_SIM` (NVT, restrained; velocities generated here) |
| 5 | NPT density equilibration (restrained, iterative, convergence-checked) |
| 6 | Relax — unrestrained NPT (optional; only if `RELAX_NS > 0`) |
| 7 | Run production MD (NPT, unrestrained) |
| 8 | Finalize outputs, create trajectory symlinks |
| 9 | Write parameters log |
| 10 | Post-analysis: PBC/strip/align + RMSD/Rg/RMSF/DSSP + clustering |

Stage folders: `build/ → em/ → heat/ → density/ → [relax/] → prod/`. By default
(`RELAX_NS=0`) restraints release at the **start of production**, so the trajectory captures
the protein relaxing away from the input pose; set `RELAX_NS > 0` to equilibrate first so
production starts pre-relaxed (e.g. for bound-state sampling). See
`scripts/simulation/MD-output-guide.md` for a full description of all output files.

### Production ensemble (`ENSEMBLE=NVT|NPT`)

T-REMD production runs **NVT** by default (constant volume). Set `ENSEMBLE=NPT` in the
submit script to pressure-couple production with a **C-rescale barostat** (at `REF_P` bar,
default 1.0) — the correct ensemble for density-sensitive observables and bound-state
sampling. Under `NPT` the per-replica `equil/` stage is also pressure-coupled, and GROMACS
automatically adds the *PV* term to the replica-exchange Metropolis criterion (`-replex` is
unchanged). The constant-temperature interpretation is unaffected either way:
`prod/rep000/remd.xtc` is still the lowest-temperature ensemble. Plain MD production is
always NPT.

**`REF_P` and `TAU_P` are not "NPT only", and they apply to all three engines.** The
`density/` stage runs the barostat in every engine and under *both* ensembles — that stage
is what determines the box size. So `REF_P` (default `1.0` bar) is the pressure the box was
equilibrated at no matter what `ENSEMBLE` is set to.

All `ENSEMBLE` controls is whether the barostat keeps running in the stages after
`density/`:

| | `density/` | `equil/` and `prod/` |
|---|---|---|
| `ENSEMBLE=NPT` | barostat on at `REF_P` | barostat on at `REF_P` |
| `ENSEMBLE=NVT` | barostat on at `REF_P` | barostat **off** — box frozen at whatever `density/` ended on |

So `REF_P` still matters under NVT: it sets the pressure the frozen box was equilibrated
at. `REF_P=10` with `ENSEMBLE=NVT` gives a constant-volume run in a box the size water
would occupy at 10 bar.

`TAU_P` (default `2.0` ps) is the C-rescale coupling time. C-rescale is thermodynamically
consistent — unlike Berendsen it reproduces the correct NPT distribution at *any*
`tau-p` — so this sets how fast the volume relaxes, not which ensemble is sampled. Keep it
at or above ~1 ps: the pressure autocorrelation time is only 0.1–0.5 ps, so a shorter
`tau-p` makes the barostat chase instantaneous virial noise, and in `density/` (where
`-DPOSRES` is on) that rescaling does spurious work against the restraints.

---

### Force-field compatibility

`FF` selects the `pdb2gmx` force field (default `amber99sb-ildn`). It accepts a short
alias from `FF_ALIASES` in `site_config.sh` — e.g. `FF="charmm36m"` — which the engine
expands to the installed dated port name and records **resolved** in `parameters.txt`,
so a finished run always names the exact release it used. Any `charmm*` name switches
the generated mdp to force-switched vdW at 1.2 nm, as CHARMM requires.

**T-REMD and plain MD accept any force field `pdb2gmx` can build.** They scale nothing,
so the topology is used exactly as written.

**REST2 does not.** It works by having `plumed partial_tempering` rewrite the topology
to scale the solute's terms by λ, and that tool only knows `[ atoms ]`, `[ atomtypes ]`,
`[ nonbond_params ]`, `[ pairtypes ]` and `[ dihedrals ]`. A force field with a solute
term outside that list is silently only *partially* scaled. CHARMM's **CMAP** backbone
cross-term is exactly such a term — measured byte-identical between λ=1.0 and λ=0.5 —
and it applies only to protein, i.e. precisely the hot region. A CHARMM REST2 run would
therefore scale charges, LJ and dihedrals while leaving the backbone conformational term
at full strength, and **neither self-check would notice**: λ=1.0 stays exact so the
`scale=1.0 → P=1.0` sanity pair passes, and exchanges still occur at normal rates so the
acceptance gate passes.

So the REST2 engine **rejects `FF=charmm*` at STEP 1** rather than produce a
plausible-looking wrong ensemble. This is a limitation of the scaling tool, not of the
GROMACS build (the 2023.5 REST2 build reads `GMXLIB` and builds CHARMM topologies fine)
and not of REST2 as a method. Note AMBER **ff19SB** also uses CMAP and would hit the
same problem; it is not guarded by name, so don't use it with REST2.

For CHARMM enhanced sampling, use T-REMD.

**Trying a force field this pipeline has not used before?**
[`docs/FORCE_FIELDS.md`](docs/FORCE_FIELDS.md) explains the rule and gives a short,
GPU-free procedure to check whether `partial_tempering` actually scales everything that
force field needs scaled — run it once before any production REST2 run.

---

## Output & Analysis

**Output model (folder-symlink).** Small, laptop-worthy dirs (`build/ em/ analysis/
logs/`, `parameters.txt`, final PDB) stay **real** in `OUTDIR`. The bulk stage dirs
(`heat/ density/ equil/ prod/` for REMD/REST2; `heat/ density/ [relax/] prod/`
for MD) are **folder symlinks into
scratch** (`SCRATCH_DIR` under `SCRATCH_ROOT`), so mdrun writes the large trajectories
straight onto scratch — the tight-quota pool is never touched and a laptop rsync of
`OUTDIR` sees only a few broken folder-links instead of thousands of file-links.
Analysis reads e.g. `prod/rep000/remd.xtc` through the symlink (no path change). On
success the run also copies the real dirs into `SCRATCH_DIR`, so the scratch archive is
self-contained. Copy anything you need long-term. Set **`SYMLINK_BULK=0`** to keep
everything **real in `OUTDIR`** with no scratch offload (e.g. when `OUTDIR` is already on
a large disk).

Post-analysis runs automatically at the end of **both** engines' jobs — PBC fix, protein
strip + backbone align, then RMSD / Rg / RMSF / DSSP and conformational clustering. T-REMD
additionally computes replica-exchange acceptance rates. The same `scripts/analysis/` tools
serve both engines (the analysis layer detects MD vs REMD automatically), and multi-chain
complexes are handled by a dedicated path. See `scripts/analysis/README.md` for the full
script reference.

**T-REMD key point:** each replica runs at a **fixed temperature** and **coordinates** are
exchanged between replicas, so:

- `prod/rep000/remd.xtc` is the lowest-temperature constant-temperature trajectory — use it directly for analysis.
- No demux step is needed.

To re-run analysis on a finished job (no simulation is re-run) — `run_analysis.sh`
auto-detects MD vs T-REMD vs REST2 from the job layout, so one command covers all three:
```bash
# Directly (login node / interactive allocation):
bash   scripts/analysis/run_analysis.sh    OUTDIR          # plain MD
bash   scripts/analysis/run_analysis.sh    OUTDIR 000      # T-REMD / REST2, slot 000
gromd-acceptance OUTDIR                                   # acceptance rates alone (target 20–30%)

# As a SLURM job (long trajectories) — copy, set OUTDIR, run:
cp example/submit_jobs/submit_analysis.sh my_analysis.sh
bash my_analysis.sh
```
Re-runs overwrite cleanly, so this is safe to repeat after editing an analysis script.
`CLUSTER_CUTOFF`, `SHELL_NM` and `N_SNAPSHOTS` are environment overrides rather than job
parameters — see `scripts/analysis/README.md`.

---

## Configuration Reference

| File | Who edits it | What it controls |
|------|-------------|-----------------|
| `site_config.sh` | Once per user/cluster | GROMACS paths, scratch root, module names, `GMXLIB`, force-field aliases |
| `my_job.sh` (copy of a `submit_*.sh`) | Per job | PDB input, force field, temperature(s), simulation length, ensemble (REMD/REST2: replicas + range) |
| `my_config.sh` (copy of `config_example.sh`) | Per job, optional | The same job parameters, in a standalone file passed as the engine's first argument |
| `#SBATCH` headers in engine script | Only if changing resource defaults | Partition, GPU type, memory, wall time |

**Full parameter reference: [`docs/PARAMETERS.md`](docs/PARAMETERS.md)** — every job
parameter, its default, its range, and which engines accept it.

---

## Directory Structure

```
gromacs_REMD/
├── site_config.sh              # Cluster-level settings (edit once)
├── pyproject.toml              # Python package: gromd_analysis + the gromd-* commands
├── docs/
│   ├── PARAMETERS.md           # Full job-parameter reference (all three engines)
│   └── FORCE_FIELDS.md         # FF selection, aliases, and the REST2 compatibility check
├── gromd_analysis/             # The Python half of the analysis layer (installable)
│   ├── layout.py                  # JobDir — parses OUTDIR into typed paths (gromd-layout)
│   ├── xvg.py                     # .xvg parser + line plot     (gromd-plot-xvg)
│   ├── dssp.py                    # secondary-structure map     (gromd-plot-dssp)
│   ├── clustering.py              # conformational clustering   (gromd-cluster)
│   ├── remd_log.py                # exchange acceptance rates   (gromd-acceptance)
│   └── chains.py                  # per-chain .ndx groups       (gromd-chain-index)
├── scripts/
│   ├── simulation/             # Simulation engines (do not edit)
│   │   ├── REMD-gromacs.sbatch    # T-REMD engine
│   │   ├── MD-gromacs.sbatch      # Plain-MD engine
│   │   ├── config_example.sh      # Job config template (copy and edit)
│   │   ├── REMD-output-guide.md   # T-REMD output file reference
│   │   └── MD-output-guide.md     # Plain-MD output file reference
│   ├── analysis/               # GROMACS-driving shell steps (see scripts/analysis/README.md)
│   │   ├── run_analysis.sh        # Whole post-analysis, re-runnable on a finished job
│   │   ├── analysis.sbatch        # Same, as a CPU-only SLURM job
│   │   └── calc_traj_*.sh, fix_PBC*.sh, multichain_*.sh   # the gmx CLI steps
│   └── installation/           # Build recipes — start at installation/README.md
│       ├── README.md              # Ordered GROMACS/PLUMED/FF/conda install guide
│       ├── install_plumed.sh      # PLUMED 2.9.4
│       ├── install_gromacs-plumed.sh          # GROMACS 2024.3+PLUMED (T-REMD/MD)
│       ├── install_gromacs-2023.5-plumed.sh   # GROMACS 2023.5+PLUMED (REST2)
│       ├── install_python_env.sh  # groMD_env + gromd_analysis
│       └── plumed.sh.template     # template for $HOME/plumed.sh
└── example/
    ├── input_pdbs/             # Example protein structures
    └── submit_jobs/
        ├── submit_REMD.sh     # T-REMD submission wrapper (reads site_config.sh)
        ├── submit_MD.sh       # Plain-MD submission wrapper (reads site_config.sh)
        └── submit_analysis.sh # Re-run analysis on an existing output dir
```
