# TODO list

- [in progress] Restructure output-file behavior so logs/metadata aren't orphaned from the trajectory data. **Design agreed — see `knowledgebase/plans/output-restructure.md`.** Approach: keep the small laptop-worthy dirs (`analysis/`, `logs/`, `em/`, `build/`, `parameters.txt`, final PDB) **real** in OUTDIR and make the bulk stage dirs (`prod/`, `equil/`, `density/`[, `heat/`, `relax/`]) **folder-level symlinks** into scratch (created before the run, so bulk writes straight to scratch → quota-safe). This gives 3–4 broken folder-links on the laptop rsync instead of thousands of file-links, keeps rerun-analysis-from-OUTDIR working via the folder symlinks, and — on success — copies the real dirs into the scratch dir so scratch is a self-contained, self-describing archive (unique jobid-based name). `build/` stays real to preserve the same-fs genion rename. Not yet implemented.

- [completed] Package the Python analysis scripts into a simple installable package. Done as a hybrid: the five GROMACS-free Python scripts became `gromd_analysis/` (flat layout, `pyproject.toml`, `gromd-*` console scripts) plus a new `layout.py`/`JobDir`; the `gmx`-driving shell scripts stayed in `scripts/analysis/`. Phase 2 (shipping the shell as package data so out-of-repo consumers get the whole pipeline from the install) is deliberately deferred.

- [pending] Replace the REMD/REST2 `equil/` stage with a per-replica `density2/`: identical to `density/` (segmented, checkpoint-chained, volume-plateau convergence via `density_converged.py`) except one instance per replica at its own `T_i` under one `mpirun -multidir`, converging when every replica plateaus — replicas end at *different* volumes by design (thermal expansion), so do not test for equal densities. NPT-only. Two open questions: whether the NVT path keeps a short fixed-length thermalization (deleting `equil/` would leave NVT with no per-replica thermalization before exchanges start), and whether the stage is restrained — see the question item below.

- [pending] **Question** — should the per-replica stage (`equil/`, or `density2/` if the item above lands) carry the same `-DPOSRES` restraints as `density/`? Restraints are currently active in `em/ heat/ density/` and NOT in `equil/ relax/ prod/` in all three engines, so REMD/REST2 release them one stage earlier than MD does. Arguments for restraining: the volume plateau test then measures box relaxation only, instead of a signal contaminated by solute conformational change at the top of the ladder; every replica enters production at the designed pose, so RMSD has one well-defined origin; and hot replicas cannot partially unfold before the first exchanges swap them down into `rep000`. Argument against: `equil/` currently absorbs the restraint-release transient before `-replex` goes live — but that relaxes in tens of ps and is largely common-mode across the ladder, so it mostly cancels in the neighbour-pair ΔE the Metropolis criterion actually uses. Settle it empirically if wanted: same system both ways, compare `rep000` RMSD over the first ns and acceptance over the first ~100 exchange attempts.

- [pending] **Question** — should analysis measure RMSD against the *designed/minimized* structure rather than frame 0 of production? `strip_and_align_trajectory.sh` currently extracts the reference with `trjconv -dump 0`, so "drift" is measured from the first production frame — which already contains every bit of drift accumulated through `heat/`, `density/` and `equil/`, and is a different structure in each replica. Using `em/em.gro` (the same anchor `POSRES_REF` uses) would make RMSD mean drift from the design, matching what these pipelines are documented to measure. The plumbing is partly there already — the script accepts a `REF_GRO` override — but note the reference must be extracted PBC-whole or every RMSD is silently inflated (see GOTCHAS).

- [pending] (optional) Switch the analysis code from MDAnalysis to mdtraj.

- [pending] **Question** — should all Python deps move into `pyproject.toml` (drop the `pip install --no-deps` habit and shrink `environment.yml` to just the interpreter + pip), so `gromd_analysis` installs into *any* conda env instead of only `groMD_env`? Confirmed pip-safe: every current dep (matplotlib, numpy, MDAnalysis, scikit-learn) and the incoming mdtraj (>=1.10) ship manylinux wheels — no conda-only package (no pyrosetta/OpenMM/GROMACS-bindings), so `--no-deps` is unnecessary for this package. *For:* portability (install anywhere), and a single source of truth — `environment.yml` and `pyproject.toml` already list deps twice and have drifted (env.yml carries `pandas`+`seaborn` that the toml lacks and the code never imports). *Against / caution:* pip's numpy is OpenBLAS vs conda's MKL, and pip upgrading numpy under a conda C-extension causes ABI skew — mitigate with loose lower bounds (`numpy>=1.23`, etc.) so pip respects an already-present conda numpy rather than rebuilding it. Not a blocker: `GROMD_ENV` in `site_config.sh` is already env-overridable, so runtime activation of any env is a separate, already-solved concern. Files it would touch: `pyproject.toml`, `scripts/installation/environment.yml`, `scripts/installation/install_python_env.sh`, `README.md`, `scripts/analysis/README.md`, plus the `--no-deps` gotcha in `CLAUDE.md` + `knowledgebase/GOTCHAS.md` (and optionally a stale comment on `site_config.sh:93`).

- [pending] add the gpu efficiency thing (cuda MMP?) to REMD. It's already implemented in the rest2 script

- [pending] **Question** - should the `density/` stage be restrained? Or should there be an unrestrained phase? relaxing the structure may have an effect on the overall volume/density and thus the pressure during the NVT phase.


- [completed] Make the EM step run with its cwd inside `em/` (or sweep `step*.pdb` afterwards). Done 2026-09-02: each engine now does `cd "${SCRATCH_DIR:-$OUTDIR}"` immediately before EM (placement is load-bearing — earlier breaks the cross-fs `solvate`/`genion` rename), and `report_crash_dumps.sh` writes a few-KB `OUTDIR/CONSTRAINT_FAILURES.txt` breadcrumb plus a job-log banner so scratch-resident dumps stay visible. When SETTLE rejects a minimization step, mdrun writes `step<N>b.pdb`/`step<N>c.pdb` with a bare relative filename, so they land in the **submit directory** — two full-system PDBs (~75 MB each on a 949k-atom system) left on the tight-quota pool. Benign output, not a crash: steepest descent rejects the step, halves `Dmax`, and continues.







## for claude
- The user will ask you to add an item to the TODO list or ask you what is on the TODO list.
- You will update the TODO list accordingly and provide the current list when requested.
- You can track the status of each item as "pending", "in progress", or "completed".
- update the conventions for this TODO list here as we decide on them.
- format - todo list items should be 1-2 sentences long, and should be clear and concise.
- an item that is an open design/science question rather than a task is written as `- [pending] **Question** — ...`, and states the arguments on each side so the decision can be made without re-deriving them.
