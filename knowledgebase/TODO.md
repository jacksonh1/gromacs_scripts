# TODO list

- [in progress] Restructure output-file behavior so logs/metadata aren't orphaned from the trajectory data. **Design agreed — see `knowledgebase/plans/output-restructure.md`.** Approach: keep the small laptop-worthy dirs (`analysis/`, `logs/`, `em/`, `build/`, `parameters.txt`, final PDB) **real** in OUTDIR and make the bulk stage dirs (`prod/`, `equil/`, `density/`[, `heat/`, `relax/`]) **folder-level symlinks** into scratch (created before the run, so bulk writes straight to scratch → quota-safe). This gives 3–4 broken folder-links on the laptop rsync instead of thousands of file-links, keeps rerun-analysis-from-OUTDIR working via the folder symlinks, and — on success — copies the real dirs into the scratch dir so scratch is a self-contained, self-describing archive (unique jobid-based name). `build/` stays real to preserve the same-fs genion rename. Not yet implemented.

- [completed] Package the Python analysis scripts into a simple installable package. Done as a hybrid: the five GROMACS-free Python scripts became `gromd_analysis/` (flat layout, `pyproject.toml`, `gromd-*` console scripts) plus a new `layout.py`/`JobDir`; the `gmx`-driving shell scripts stayed in `scripts/analysis/`. Phase 2 (shipping the shell as package data so out-of-repo consumers get the whole pipeline from the install) is deliberately deferred.

- [pending] (optional) Switch the analysis code from MDAnalysis to mdtraj.

- [pending] add the gpu efficiency thing (cuda MMP?) to REMD. It's already implemented in the rest2 script

- [pending] Make the EM step run with its cwd inside `em/` (or sweep `step*.pdb` afterwards). When SETTLE rejects a minimization step, mdrun writes `step<N>b.pdb`/`step<N>c.pdb` with a bare relative filename, so they land in the **submit directory** — two full-system PDBs (~75 MB each on a 949k-atom system) left on the tight-quota pool. Benign output, not a crash: steepest descent rejects the step, halves `Dmax`, and continues.







## for claude
- The user will ask you to add an item to the TODO list or ask you what is on the TODO list.
- You will update the TODO list accordingly and provide the current list when requested.
- You can track the status of each item as "pending", "in progress", or "completed".
- update the conventions for this TODO list here as we decide on them.
- format - todo list items should be 1-2 sentences long, and should be clear and concise.
