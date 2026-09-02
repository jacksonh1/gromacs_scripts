# Installation — GROMACS, PLUMED, force fields, Python env

Everything the pipeline needs, in the order it has to be built. Read this once per
cluster; after it, `site_config.sh` is the only file you touch.

There are **two independent halves**:

| Half | Needed for | Where it goes |
|------|-----------|---------------|
| GROMACS + PLUMED builds | running any simulation | `$HOME/opt/{gromacs,plumed}/…` |
| `groMD_env` conda env | post-analysis / plotting | conda envs |

If someone else on your cluster already built GROMACS, skip to
[Step 7](#step-7--configure-site_configsh) and just point `site_config.sh` at it.

---

## What you end up with

```
$HOME/opt/
├── plumed/2.9.4/                  # PLUMED (kernel + CLI)
├── gromacs/2024.3-plumed/         # main build — T-REMD and plain MD  (GMXRC)
├── gromacs/2023.5-plumed/         # REST2 only — working PLUMED -hrex (REST2_GMXRC)
├── gromacs/ff/<charmm-port>.ff/   # optional: CHARMM force field      (GMXLIB)
└── hcoll_compat/                  # 2 staged .so files                (HCOLL_COMPAT_DIR)
$HOME/plumed.sh                    # env script sourced by the REST2 engine (PLUMED_SH)
```

**You only need the 2023.5 build if you intend to run REST2.** T-REMD and plain MD
run entirely on `2024.3-plumed`.

### Why two GROMACS builds

The PLUMED `-hrex` (Hamiltonian replica exchange) patch was never ported to the
GROMACS-2024 series: it patches, runs, and reports success while producing
`dE_term = 0` and **zero exchanges** — a silently fake REST2. GROMACS 2023.5 is the
newest release with working hrex. T-REMD and MD don't use `-hrex`, so they stay on
2024.3. See the "-hrex is SILENTLY BROKEN" entry in `knowledgebase/GOTCHAS.md`.

### Why the `-plumed` build is the default even for T-REMD

`site_config.sh` points `GMXRC` at `2024.3-plumed` because it is the build that
exists and is validated here. A plain (unpatched) 2024.3 works identically for
T-REMD/MD — `install_gromacs.sh` builds one — but nothing in the pipeline needs it.
Treat `install_gromacs.sh` as a spare, not part of the path below.

---

## Before you start

- **Node choice matters.** GROMACS bakes the build host's CPU SIMD level and
  compiler `-march` into the binary. Build on the *oldest* microarch you will run
  on, or force the target explicitly (Step 4 does the former, Step 5 the latter).
  On this cluster the run nodes `node[3619-3620]` are Skylake-AVX512.
- **All build scripts are `sbatch` scripts** with `#SBATCH -o ./logs/…`, so submit
  them **from the repo root** (which has `logs/`) or any dir containing `logs/`.
- **Sources are not downloaded by the scripts.** Each script `cd`s into an
  already-unpacked source tree. Download and unpack them first (Step 1 / Step 3).
- Pick a directory to hold source trees. On this cluster that is
  `/orcd/pool/004/$USER/installations/`; anywhere with ~20 GB free works.
- Budget: PLUMED ~20 min, each GROMACS ~1–2 h with `make check`.

---

## Step 1 — Build PLUMED 2.9.4

```bash
cd /orcd/pool/004/$USER/installations
wget https://github.com/plumed/plumed2/releases/download/v2.9.4/plumed-2.9.4.tgz
tar xzf plumed-2.9.4.tgz          # -> ./plumed-2.9.4/
```

`install_plumed.sh` expects `./plumed-2.9.4/` relative to **its own submit
directory**, so submit it from the directory holding the source tree — and that
directory also needs a `logs/`:

```bash
mkdir -p logs
sbatch /path/to/gromacs_REMD/scripts/installation/install_plumed.sh
```

It configures with `--prefix=$HOME/opt/plumed/2.9.4`, then `make -j8`,
`make regtest`, `make install`. Version pinning is deliberate: PLUMED 2.9.x is the
series verified against GROMACS 2023.5 hrex.

Check: `$HOME/opt/plumed/2.9.4/bin/plumed info --version` → `2.9.4`.

## Step 2 — Write `$HOME/plumed.sh`

The GROMACS build scripts and the REST2 engine both `source ~/plumed.sh` to put
PLUMED on `PATH` and expose the kernel. Copy the shipped template and edit the two
paths at the top:

```bash
cp scripts/installation/plumed.sh.template ~/plumed.sh
$EDITOR ~/plumed.sh        # set PLUMED_PREFIX and PLUMED_SRC
```

`site_config.sh` finds it via `PLUMED_SH` (default `$HOME/plumed.sh`), so you can
keep it elsewhere if you prefer.

## Step 3 — Download the GROMACS sources

```bash
cd /orcd/pool/004/$USER/installations
wget https://ftp.gromacs.org/gromacs/gromacs-2024.3.tar.gz && tar xzf gromacs-2024.3.tar.gz
wget https://ftp.gromacs.org/gromacs/gromacs-2023.5.tar.gz && tar xzf gromacs-2023.5.tar.gz   # REST2 only
```

Keep each tree **clean and unpatched** until the build script patches it. If you
ever need to start over, re-extract rather than trying to unpatch.

## Step 4 — Build GROMACS 2024.3 + PLUMED (T-REMD / MD)

Submit from the directory containing `gromacs-2024.3/` (plus a `logs/`):

```bash
sbatch /path/to/gromacs_REMD/scripts/installation/install_gromacs-plumed.sh
```

What it does: loads `openmpi/5.0.8` + `cuda/12.9.1`, sources `~/plumed.sh`, runs
`plumed patch -p -e gromacs-2024.3` **inside the job**, then cmake with
`-DGMX_MPI=on -DGMX_THREAD_MPI=OFF -DGMX_GPU=CUDA -DGMX_BUILD_OWN_FFTW=ON`,
`make -j8`, `make check`, `make install` into `$HOME/opt/gromacs/2024.3-plumed`.

Two things to know:

- The in-job `plumed patch` only works on a node where the `plumed` **CLI** can run.
  On `pi_keating` nodes it can; on some `mit_normal_gpu` nodes it dies on a missing
  `libhcoll.so.1` (Step 6). If the job fails at the patch step, patch on the login
  node first and delete the `plumed patch` line, exactly as the 2023.5 script does.
- The script has **no `set -e`**, so a failing `make` still falls through to
  `make install`. Read the log and confirm `Installing …` lines before trusting it.

Check: `$HOME/opt/gromacs/2024.3-plumed/bin/gmx_mpi --version` — note there is
**only `gmx_mpi`**, no plain `gmx`, throughout this project.

## Step 5 — Build GROMACS 2023.5 + PLUMED (REST2 only)

Skip unless you will run REST2.

This build differs from Step 4 in four ways: cuda **12.4.0** (2023.5 predates 12.9),
`cmake/3.24.3` (cmake 4.x rejects 2023.5's bundled googletest), forced
`-DGMX_SIMD=AVX_512 -march=skylake-avx512` so it can build on any free node yet still
run on the Skylake run nodes, and **pre-patching on the login node**.

Patch first, on the login node:

```bash
export SRC=/orcd/pool/004/$USER/installations/gromacs-2023.5
source ~/plumed.sh
cd "$SRC" && rm -rf build && plumed patch -p -e gromacs-2023.5
```

Then submit (the script reads `SRC`, defaulting to the path above with `jhalpin`
hardcoded — export your own):

```bash
cd /path/to/gromacs_REMD
SRC=$SRC sbatch --export=ALL scripts/installation/install_gromacs-2023.5-plumed.sh
```

The script refuses to build if it doesn't find `.preplumed` backups in the source,
and it tolerates `make check` failures (regression tests can fail on unrelated GPU
precision checks) because the authoritative test is the hrex acceptance check below.

**Verify hrex actually exchanges** before trusting any REST2 result — a build with
broken hrex looks identical to a working one until you count swaps. Run a short
REST2 job and check the acceptance rate is nonzero (a `scale=1.0` sanity pair should
give P ≈ 1.0). The REST2 engine has a self-check for this; don't disable it.

## Step 6 — Stage the hcoll compat libraries

`gmx_mpi` and the `plumed` CLI carry a baked-in `NEEDED` on `libhcoll.so.1` and
`libocoms.so.0` from `/opt/mellanox/hcoll/lib`. Those libs are **inert at runtime**
(the loaded OpenMPI routes around them) but the loader still must resolve them at
exec, and `/opt/mellanox` is node-local — missing on some GPU nodes, where the job
dies instantly with `error while loading shared libraries: libhcoll.so.1`.

From a node that has them (the login node does):

```bash
mkdir -p ~/opt/hcoll_compat
cp -a /opt/mellanox/hcoll/lib/libhcoll.so.1* \
      /opt/mellanox/hcoll/lib/libocoms.so.0* ~/opt/hcoll_compat/
```

The engines prepend `HCOLL_COMPAT_DIR` to `LD_LIBRARY_PATH`; it's a harmless no-op
on nodes that already have the libs. If your cluster's GROMACS doesn't link them
(`ldd $(which gmx_mpi) | grep hcoll` finds nothing), skip this step entirely.

Related: the engines set `#SBATCH -C rocky8` so the binary can't land on a node with
a different glibc than the one it was built against.

## Step 7 — Install CHARMM36m (optional)

Only needed if you want to run CHARMM36m. GROMACS bundles `charmm27.ff`, not 36m.

Download the GROMACS port of CHARMM36m from the MacKerell lab page
(<http://mackerell.umaryland.edu/charmm_ff.shtml>, "CHARMM36 force field in GROMACS
format"). Take the **force-switch** port — **not** the `ljpme` one, which needs
`vdwtype = PME` and a different mdp path than this pipeline generates.

```bash
mkdir -p ~/opt/gromacs/ff
tar xzf charmm36-*.ff.tgz -C ~/opt/gromacs/ff     # -> ~/opt/gromacs/ff/charmm36-feb2026_cgenff-5.0.ff
```

**Do not rename the extracted directory.** These ports ship with the release date in
the name (`charmm36-jul2022.ff`, `charmm36-feb2026_cgenff-5.0.ff`, …) and you pass
that name, minus `.ff`, to `pdb2gmx -ff`. A generic name like `charmm36m.ff` reads
better in a submit script, but then `parameters.txt` records `FF=charmm36m` forever
— and dropping a newer port into that same directory later silently changes the
force field of every future run while the job record looks unchanged. The mdp gate is
a `charmm*` glob, so the long dated name needs no code change.

Because that name is a mouthful, `site_config.sh` carries an alias table
(`FF_ALIASES`) mapping a short name to the installed directory. The engine resolves
the alias at STEP 1 and writes the **resolved** name into `parameters.txt`, so the
convenience costs no provenance. Add your port to the table:

```bash
declare -A FF_ALIASES=(
  [charmm36m]="charmm36-feb2026_cgenff-5.0"
)
```

Then a job can use either name:

```bash
FF="charmm36m"    # in your submit script; the real directory name also works
WATER="tip3p"
```

An alias pointing at a force field that isn't installed fails at STEP 1 with a clear
message rather than deep inside `pdb2gmx`.

> **Installing a force field other than CHARMM36m?** Two things to check before using
> it, both covered in [`../../docs/FORCE_FIELDS.md`](../../docs/FORCE_FIELDS.md):
> its published **nonbonded settings** (the mdp generator only knows the AMBER
> plain-cutoff and CHARMM force-switch cases, and silently gives anything else the
> AMBER treatment), and — if you intend to run **REST2** — whether
> `plumed partial_tempering` actually scales every solute term it uses. A force field
> it does not fully scale produces a healthy-looking run of the wrong ensemble.

`site_config.sh` exports `GMXLIB=$HOME/opt/gromacs/ff`, which GROMACS searches *in
addition to* each build's bundled `top/`, so AMBER keeps working. Keep
`WATER=tip3p`: inside a CHARMM FF dir that resolves to the CHARMM-modified TIP3P
automatically. Never cross a force field with another one's water.

## Step 8 — Configure `site_config.sh`

Edit the repo-root `site_config.sh`. Every value is `${VAR:-default}`, so an
environment variable can override any of them per job.

| Variable | Set it to |
|---|---|
| `GMXRC` | `$HOME/opt/gromacs/2024.3-plumed/bin/GMXRC` (Step 4) |
| `REST2_GMXRC` | `$HOME/opt/gromacs/2023.5-plumed/bin/GMXRC` (Step 5) |
| `GMXLIB` | dir holding the CHARMM `.ff` port (Step 7) |
| `PLUMED_SH` | your `plumed.sh` (Step 2) |
| `HCOLL_COMPAT_DIR` | staged libs (Step 6) |
| `SCRATCH_ROOT` | large-disk root for trajectories, ~100 GB per job |
| `CUDA_MODULE`, `OPENMPI_MODULE` | modules matching the Step 4 build |
| `REST2_CUDA_MODULE` | module matching the Step 5 build (`cuda/12.4.0`) |
| `CONDA_MODULE`, `GROMD_ENV` | conda module and analysis env name |

The CUDA/MPI modules must match what each build was compiled against — that is why
REST2 has its own `REST2_CUDA_MODULE`.

## Step 9 — Create the analysis conda env

```bash
bash scripts/installation/install_python_env.sh
```

Builds `groMD_env` from `environment.yml` and installs `gromd_analysis` editable
into it, which is what puts the `gromd-*` commands on `PATH`. The pipeline activates
this env only for the post-analysis steps — the simulation steps use the system
`python3` (stdlib only). Details and the no-install fallback are in the repo README.

Check: `conda activate groMD_env && pytest` from the repo root.

---

## Verifying the whole install

```bash
source site_config.sh
source "$GMXRC"
gmx_mpi --version | head -20            # GROMACS + CUDA + MPI + SIMD line
ls "$GMXLIB"/*.ff                       # force fields installed in Step 7

# Step 7 end-to-end: build a topology with the CHARMM port (in a throwaway dir)
gmx_mpi pdb2gmx -f example/input_pdbs/helix_fusion.pdb -o t.gro -p t.top -i p.itp \
    -ff charmm36-feb2026_cgenff-5.0 -water tip3p -ignh
# -> "Using the Charmm36-feb2026_cgenff-5.0 force field in directory $GMXLIB/..."
```

Note there is **no `-ff list`** in GROMACS 2024 — it reads `list` as a force-field
name and aborts. To browse what's available, run `pdb2gmx` with no `-ff` and read the
interactive menu, or just `ls "$GMXLIB"/*.ff` plus the build's own
`share/gromacs/top/`.

Then run a real short job end to end — a 2 ns T-REMD on
`example/input_pdbs/` — which exercises build, EM, equilibration, production,
and the whole analysis layer. That is the only check that covers everything.

---

## Porting to another cluster

The build scripts are written for the Keating lab's MIT SuperCloud/ORCD partitions.
Every site-specific line is in the `#SBATCH` header or the module block:

| What | Where | Change to |
|---|---|---|
| `-p pi_keating` / `-p mit_normal_gpu` | all four scripts | your partition |
| `--nodelist=node[3619-3620]` | `install_plumed.sh`, `install_gromacs*.sh` | your oldest-microarch build node, or drop it and pin SIMD instead |
| `--constraint="rocky8"` | all | your OS constraint, or remove |
| `--gres=gpu:1` / `gpu:l40s:1` | GROMACS scripts | a GPU matching your run target (used by `make check` and CUDA arch selection) |
| `module load openmpi/5.0.8`, `cuda/12.9.1`, `cuda/12.4.0`, `cmake/3.24.3-x86_64` | all | your module names |
| `-DGMX_SIMD=AVX_512`, `-march=skylake-avx512` | `install_gromacs-2023.5-plumed.sh` | your oldest run-node microarch, or delete both to auto-detect |
| `SRC=/orcd/pool/004/jhalpin/installations/gromacs-2023.5` | `install_gromacs-2023.5-plumed.sh` | your source tree (or `export SRC=` at submit) |
| `$HOME/opt/...` prefixes | all | wherever you install |
| `export LD_LIBRARY_PATH=` (clearing it) | all | keep — it prevents a stale conda/module `LD_LIBRARY_PATH` from leaking into the build |
| hcoll staging (Step 6) | — | skip if `ldd` shows no hcoll dependency |

The **rules** that survive the port, whatever the cluster:

1. Build on the oldest CPU microarch you will run on, or force `-march`/`GMX_SIMD`.
2. Build PLUMED before GROMACS; patch the GROMACS source before cmake.
3. Match each GROMACS build's CUDA/MPI modules at run time.
4. Never use a GROMACS-2024 build for REST2.

---

## Script inventory

| Script | Purpose | Notes |
|---|---|---|
| `install_plumed.sh` | PLUMED 2.9.4 | needed for both GROMACS builds |
| `install_gromacs-plumed.sh` | GROMACS 2024.3 + PLUMED → `GMXRC` | **the main build**; T-REMD + MD |
| `install_gromacs-2023.5-plumed.sh` | GROMACS 2023.5 + PLUMED → `REST2_GMXRC` | REST2 only; requires login-node pre-patch |
| `install_gromacs.sh` | plain GROMACS 2024.3, no PLUMED | spare; not used by the pipeline |
| `install_python_env.sh` | `groMD_env` + `gromd_analysis` | login node, no SLURM job |
| `environment.yml` | conda env spec | — |
| `plumed.sh.template` | template for `$HOME/plumed.sh` | edit the two paths at the top |

Deeper background for every "why" above lives in `knowledgebase/GOTCHAS.md`.
