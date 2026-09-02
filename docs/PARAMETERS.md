# Job parameters

Every knob a **job** can set, for all three engines. Cluster-level settings live in
`site_config.sh` instead and are documented in
[`scripts/installation/README.md`](../scripts/installation/README.md); analysis knobs
are in [`scripts/analysis/README.md`](../scripts/analysis/README.md).

The authoritative list is `scripts/simulation/validate_params.py` (names, types,
ranges) plus the `${VAR:-default}` reads at STEP 1 of each engine. This page is the
readable version of both — if they disagree, the code wins, and this page is a bug.

---

## How to set parameters

Two ways, and they can be mixed. Both go through the same validation.

### 1. In the submit script (the common case)

Copy a template from `example/submit_jobs/`, edit the assignments at the top, run it.
The script passes them to `sbatch --export`:

```bash
cp example/submit_jobs/submit_REMD.sh my_job.sh
$EDITOR my_job.sh        # set PDB_IN, REPLICAS, T_MAX, TOTAL_NS, FF, …
bash my_job.sh
```

### 2. In a separate config file

All three engines take a config file as their **first argument** and `source` it. Use
this when you want the parameters versioned separately from the submission, or want to
run the same config through several engines:

```bash
cp scripts/simulation/config_example.sh my_config.sh
$EDITOR my_config.sh
bash example/submit_jobs/submit_REMD_with_config.sh my_config.sh
# extra sbatch flags may follow the config path:
bash example/submit_jobs/submit_REMD_with_config.sh my_config.sh -p other_partition
```

Or straight to the engine:

```bash
sbatch -n 48 scripts/simulation/REMD-gromacs.sbatch my_config.sh
```

`submit_REMD_with_config.sh` is the only shipped wrapper, but `MD-gromacs.sbatch` and
`REST2-gromacs.sbatch` accept a config argument identically — copy the wrapper and
change the engine path.

**Precedence:** the config file is sourced *after* the environment is inherited, so a
value in the config file **overrides** the same name passed via `--export`. Anything set
in neither place falls back to the engine default in the tables below.

### Validation (STEP 1, before any GROMACS work)

- **Unknown keys** — only checked for a **config file**: a line like `REPLICASS=48` is
  rejected by name before the file is sourced. Under `--export=ALL` the job inherits the
  whole shell environment, so there is no manifest to check a typo against; a misspelled
  variable set in a submit script is simply never read.
- **Values** — checked for both paths: `PDB_IN` must exist, numbers must parse, ranges
  must be sane (`T_MAX > T_MIN`, `REPLICAS >= 2`, `ENSEMBLE` in `NVT|NPT`, …).

Failures print `[ERROR] …` and stop the job immediately. Adding a parameter means
updating **three** places in sync — the engine's `${VAR:-default}` read, the `export`
list, and `validate_params.py`'s allowlist + rules.

---

## Shared parameters (all three engines)

| Parameter | Default | Meaning |
|---|---|---|
| `PDB_IN` | *(required)* | Input structure. Must exist. Always a **folded/designed pose** — it is the reference all analysis is measured against. |
| `OUTBASE` | basename of `PDB_IN` | Short name used in filenames. |
| `OUTDIR` | `remd_`/`md_`/`rest2_` + `$OUTBASE` | Output directory. **Must be on the same filesystem as the submit directory** (`solvate`/`genion` do a cross-fs-fatal `rename()`). |
| `FF` | `amber99sb-ildn` | `pdb2gmx` force field. Accepts an alias from `FF_ALIASES` in `site_config.sh` (e.g. `charmm36m`); the engine resolves it and records the real dated name in `parameters.txt`. See [Force fields](#force-fields). |
| `WATER` | `tip3p` | Water model, resolved **inside** the FF directory, so it always matches the force field. Do not cross them. |
| `BOX_SHAPE` | `dodecahedron` | `editconf -bt`. Dodecahedron ≈ truncated octahedron, ~29 % fewer waters than cubic. |
| `BOX_BUFFER` | `1.0` | nm from solute to box edge (1.0 nm = 10 Å). |
| `NEUTRALIZE` | `1` | `0`/`1`. Add counter-ions to neutralize net charge. |
| `SALT_MOLAR` | `0.15` | Additional NaCl, mol/L. `>= 0`. |
| `DT_PS` | `0.002` | Integration timestep, ps. |
| `CUTOFF_NM` | `1.0` | Non-bonded cutoff, nm. **Forced to `1.2` when `FF=charmm*`** (with an `[INFO]` notice). |
| `GAMMA_LN` | `2.0` | Thermostat coupling `tau_t`, ps⁻¹. |
| `TOTAL_NS` | `20` (REMD/REST2), `100` (MD) | Production length per replica / per trajectory, ns. |
| `SEED` | `-1` | `-1` = fresh random seed per run. `>= 0` pins `gen-seed`/`ld-seed`/`mdrun -reseed` (per-replica `SEED+i`). **Does not give bit-for-bit GPU trajectories** — it pins the setup and RNG streams only. |
| `DENSITY_SEG_STEPS` | `10000` | Steps per NPT density-equilibration segment. |
| `DENSITY_MIN_SEG` | `8` | Minimum segments before the convergence check. Must be `<= DENSITY_MAX_SEG`. |
| `DENSITY_MAX_SEG` | `20` | Maximum segments. |
| `DENSITY_TOL_REL` | `0.005` | Relative volume-change tolerance for convergence. |
| `SYMLINK_BULK` | `1` | `1` = bulk stage dirs are folder symlinks into scratch (quota-safe). `0` = everything real in `OUTDIR`, no scratch offload. |
| `SCRATCH_ROOT` | `/tmp` (overridden in `site_config.sh`) | Root for the per-job scratch directory. ~100 GB per job. |
| `SCRATCH_DIR` | derived from `SCRATCH_ROOT` | Override the exact scratch path. |
| `PRESERVE_SCRATCH_FROM` | `prod` | `prod`\|`density`\|`always`\|`never` — from which stage on scratch is kept. Stage **names**, so the meaning is the same across engines. |
| `GMX` | `gmx_mpi` | GROMACS binary. This build ships only `gmx_mpi` — there is no plain `gmx`. |

## T-REMD only (`REMD-gromacs.sbatch`)

| Parameter | Default | Meaning |
|---|---|---|
| `REPLICAS` | SLURM `-n`, else `48` | Number of temperature slots. `>= 2`. **SLURM's `-n` wins**; a config `REPLICAS` that disagrees is a hard error telling you which to change. |
| `T_MIN` | `300` | Lowest replica temperature, K. `rep000` is this ensemble. |
| `T_MAX` | `400` | Highest replica temperature, K. Must be `> T_MIN`. |
| `TEMPS_LIST` | *(unset)* | Explicit ladder, e.g. `"300.0,305.2,310.5,…"`, overriding the computed geometric one. |
| `REPLEX_PS` | `1.0` | Exchange-attempt interval, ps. **Never set below 1.0** — see the warning below. |
| `EQUIL_NS` | `0.2` | Per-replica equilibration before production, ns. `0` allowed. |
| `ENSEMBLE` | `NVT` | `NVT` (constant volume) or `NPT` (C-rescale barostat). Under NPT, GROMACS adds the *PV* term to the exchange criterion automatically. |
| `REF_P` | `1.0` | Reference pressure, bar (NPT only). |
| `TAU_P` | `1.0` | Barostat coupling time, ps (NPT only). |
| `NTOMP_SERIAL` | `min(node cores, 8)` | OpenMP threads for the single-rank setup steps (EM, density). |
| `FORCE` | `0` | `0`/`1`. Overwrite an existing `OUTDIR`. |

> **`REPLEX_PS` below 1.0 deadlocks GPU-resident REMD.** At 0.5 ps the job hangs
> mid-run with clean physics (CPUs 100 %, GPUs 0 %) — an MPI collective deadlock at
> exchange, not a blow-up. It is a threshold, not a gradient: 1.0 never hangs. The
> default is 1.0 for this reason; sub-ps exchange buys no extra sampling anyway. If you
> genuinely need it, add `-update cpu` to the production `mdrun` (~10–25 % slower).

## Plain MD only (`MD-gromacs.sbatch`)

| Parameter | Default | Meaning |
|---|---|---|
| `T_SIM` | `300` | Production temperature, K. |
| `TRAJ_PS` | `10` | Trajectory write interval, ps. |
| `HEAT_NS` | `0.2` | NVT thermalization (restrained; velocities generated here), ns. `0` allowed. |
| `RELAX_NS` | `0` | Optional unrestrained NPT relaxation before production, ns. `0` (default) releases restraints **at the start of production**, so the trajectory captures the protein relaxing away from the input pose. Set `> 0` to start production pre-relaxed — e.g. for bound-state sampling. |
| `NTOMP` | `$SLURM_CPUS_PER_TASK`, else `8` | OpenMP threads. `>= 1`. |

Production is always NPT.

## REST2 only (`REST2-gromacs.sbatch`)

REST2 takes the **same parameter set as T-REMD**, with the same types and ranges. Two
differences in meaning, and one restriction:

| Parameter | Default | Difference from T-REMD |
|---|---|---|
| `T_MIN` | `300` | The single **physical** temperature every replica runs at. |
| `T_MAX` | `450` | The maximum **effective solute** temperature. The ladder scales the solute's force field, not the thermostat. |
| `REPLICAS` | SLURM `-n`, else `24` | Fewer replicas than T-REMD reach the same range, because only the solute is heated. |
| `REPLEX_PS` | `1.0` | Already at the safe value. |
| `FF` | `amber99sb-ildn` | **CHARMM force fields are rejected** — see below. |

### Force fields are not all REST2-compatible

REST2 works by having `plumed partial_tempering` rewrite the processed topology to
scale the solute's terms by λ. It knows how to scale `[ atoms ]`, `[ atomtypes ]`,
`[ nonbond_params ]`, `[ pairtypes ]` and `[ dihedrals ]` — and **nothing else**. Any
force field whose solute energy includes a term outside that list is silently only
*partially* scaled.

CHARMM is the case that bites here: its **CMAP** backbone φ/ψ cross-term
(`[ cmap ]` / `[ cmaptypes ]`) is passed through untouched — measured byte-identical
between λ=1.0 and λ=0.5. CMAP applies only to protein, i.e. exactly the hot region, so
a CHARMM REST2 run would scale charges, LJ and dihedrals while leaving the backbone
conformational term at full strength. Neither self-check notices: λ=1.0 is still exact
so the `scale=1.0 → P=1.0` sanity pair passes, and exchanges still occur at normal
rates so the acceptance gate passes.

**So the engine refuses `FF=charmm*` at STEP 1.** This is a limitation of the *scaling
tool*, not of the GROMACS build or of REST2 as a method — the REST2 GROMACS 2023.5
build reads `GMXLIB` and builds CHARMM topologies fine.

- **AMBER** (`amber99sb-ildn`, `amber14sb`) — no CMAP. Supported.
- **CHARMM** (`charmm36m`, any `charmm*`) — CMAP. Rejected.
- **AMBER ff19SB** — also uses CMAP. It would hit the same problem and is not currently
  guarded by name; do not use it with REST2.
- Lifting the restriction means teaching `partial_tempering` to scale the CMAP grids,
  then re-validating.

For CHARMM enhanced sampling, use **T-REMD** — it scales nothing, so every force field
is handled exactly as `pdb2gmx` wrote it.

**Before running REST2 with any force field not listed above, run the compatibility
check in [`FORCE_FIELDS.md`](FORCE_FIELDS.md).** It takes a couple of minutes, needs no
GPU, and is the only way to find this class of problem — a partially-scaled run looks
healthy in every log the pipeline produces.

---

## Force fields

Set `FF` to either an alias or the real force-field directory name minus `.ff`.
Aliases live in `FF_ALIASES` in `site_config.sh`:

```bash
FF="charmm36m"                       # alias
FF="charmm36-feb2026_cgenff-5.0"     # the real directory name — equivalent
```

The engine resolves the alias at STEP 1, logs the expansion, and writes the **resolved**
name into `parameters.txt`, so a finished run always records the exact release it used.
Force fields bundled with GROMACS (`amber99sb-ildn`, `amber14sb`, …) are used by name
directly. Installing a new one is Step 7 of
[`scripts/installation/README.md`](../scripts/installation/README.md).

Full detail — aliases, water resolution, adding a force field, and the REST2
compatibility check — is in [`FORCE_FIELDS.md`](FORCE_FIELDS.md).

Any `FF` starting with `charmm` switches the generated mdp to force-switched vdW
(`vdw-modifier=force-switch`, `rvdw-switch=1.0`, `DispCorr=no`) and forces
`CUTOFF_NM=1.2`, because CHARMM was parameterized that way. The AMBER path is
unchanged and byte-identical to before that gating existed.

`WATER` is resolved inside the FF directory, so `tip3p` means standard TIP3P under
AMBER and CHARMM-modified TIP3P under CHARMM. Never pair a force field with another
one's water.

---

## What is *not* a job parameter

- **Cluster settings** — GROMACS paths, modules, scratch root, `GMXLIB`, force-field
  aliases: `site_config.sh`, edited once. A job config *may* override them (they are on
  the allowlist so they aren't flagged as typos), but normally shouldn't.
- **SLURM resources** — partition, GPU type, memory, wall time: the `#SBATCH` headers in
  the engine, or flags passed through the `_with_config` wrapper.
- **Analysis knobs** — `CLUSTER_CUTOFF`, `SHELL_NM`, `N_SNAPSHOTS` are environment
  overrides on the analysis scripts, deliberately kept out of the validated job
  parameter set. See `scripts/analysis/README.md`.
