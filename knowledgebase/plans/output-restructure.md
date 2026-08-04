# Plan: output-file restructure (folder-symlink model)

**Status:** IMPLEMENTED (2026-07-29). Supersedes TODO item 1 ("Restructure output-file
behavior so logs/metadata aren't orphaned from the trajectory data"). Implemented across all
three engines (REMD/MD/REST2) + the restart engine, with a shared
`scripts/simulation/finalize_outputs.sh`. **Added beyond the original plan:** a `SYMLINK_BULK`
toggle (default 1 = folder-symlink model below; 0 = everything real in OUTDIR, no scratch
offload) wired through the validator, both `config_example.sh`, and all engines. The
folder-symlink contract is documented in `CLAUDE.md` ("Output model" subsection) and
`knowledgebase/GOTCHAS.md`; the three output-guide docs, both READMEs, and
`REMD_log_reference.md` were updated (no `trajectories/` dir). Verified: `bash -n` on all
engines, an isolated `finalize_outputs.sh` smoke test, and `run_analysis.sh` end-to-end on the
example T-REMD dir through the new `prod/rep000/remd.xtc` path.

## Problem

Today the run is split across two locations and neither is self-contained:
- **OUTDIR** (under the repo pool `/orcd/pool`, **tight quota**, rsync'd to the user's laptop)
  holds all small files real, plus `trajectories/` full of **per-file symlinks** into scratch.
- **SCRATCH** (`/orcd/data/keating/001/$USER/MD`, a **large, durable, NOT auto-wiped** drive)
  holds only the `.xtc` trajectories (flat, `remd_rep000.xtc`) + a copy of `parameters.txt`.

Pain points this causes:
1. The laptop rsync drags **thousands of broken per-file symlinks** (`trajectories/` + the
   in-rep `prod/rep*/remd.xtc` links). Symlinks resolve on the cluster but are broken on the laptop.
2. The bulky replicated per-rep intermediates (`.tpr/.cpt/.edr/.gro` × N reps × equil+prod) sit
   **real on the tight-quota pool**.
3. Neither OUTDIR nor SCRATCH is a self-contained record: SCRATCH lacks logs/metadata/analysis;
   OUTDIR's trajectories are just links.

**Key infra facts** (see memory `project_filesystems`): the pool is durable but tight; the data
drive ("scratch") is a big drive, **not** fast and **not** periodically wiped. So the reason to
offload is **quota + laptop-sync hygiene**, not speed.

## The model: folder-level symlinks

Keep the small, laptop-worthy dirs **real** in OUTDIR; make the **bulk stage dirs symlinks**
into scratch. OUTDIR and SCRATCH then share the same tree; the only difference is which top-level
stage dirs are real vs symlinked.

| OUTDIR entry | Real / symlink | Rationale |
|---|---|---|
| `analysis/`, `logs/`, `em/`, `build/`, `parameters.txt`, `<base>_final*.pdb` | **real** | small, openable on the laptop; `build/` real also keeps the genion rename same-fs (see gotcha) |
| `prod/`, `equil/`, `density/` (REMD) — plus `heat/`, `relax/` (MD) | **symlink → `SCRATCH_DIR/<stage>`** | all bulk; created *before* the run so mdrun writes straight to scratch |

Result:
- **Rerun analysis from OUTDIR** ✓ — on the cluster `run_analysis.sh OUTDIR` reads
  `prod/rep000/remd.{tpr,xtc}` through the `prod/` folder symlink. No path change needed there.
- **Laptop sync** ✓ — exactly **3–4 broken folder-links** (`prod`, `equil`, `density`[, `heat`,
  `relax`]) instead of thousands of file-links. Harmless to leave, or one `rsync --exclude` each.
- **Bulk never touches the pool** ✓ — the folder symlinks pre-exist, so the bulk stages write
  directly onto scratch during the run. Quota-safe.
- **SCRATCH self-contained** ✓ — the bulk stages are already real there; on success we **copy**
  the real OUTDIR dirs (`analysis/`, `logs/`, `em/`, `build/`, `parameters.txt`, final PDB) into
  `SCRATCH_DIR` so the archive stands alone and is self-describing (each scratch dir carries its
  own `parameters.txt`). SCRATCH subdir name stays the **unique** `${SLURM_JOB_ID}_${TS}` (a
  browsing script identifies runs by reading `parameters.txt`, not by a descriptive dir name).

## Locked decisions

- Symlink at **folder** granularity, not per-file. Only the bulk stage dirs.
- `build/` and `em/` stay **real in OUTDIR** — small, and keeping `build/` on the pool preserves
  the same-filesystem `solvate`/`genion` topology rename (documented gotcha; no `cd` fix needed).
- SCRATCH subdir keeps a **unique, non-descriptive** name; self-description comes from the copied
  `parameters.txt` (and the copied `analysis/`).
- The current `OUTDIR/trajectories/` collection folder is **dropped** — the stage folder symlinks
  make it redundant, and it was a second source of broken laptop links.
- One shared finalize helper (`scripts/simulation/finalize_outputs.sh`) called by both engines,
  rather than duplicated logic — isolates the complexity (per repo convention).

## Concrete changes

### Both engines — before the run (replace the tree-creation line)
- REMD line 187 / MD line 169 currently `mkdir -p "${OUTDIR}"/{...,trajectories} ./logs`.
  New: create the **real** dirs in OUTDIR (`build em logs analysis` [+ `em`]) and, for each **bulk**
  stage, `mkdir -p "${SCRATCH_DIR}/<stage>"` then `ln -s "${SCRATCH_DIR}/<stage>" "${OUTDIR}/<stage>"`.
  Drop `trajectories/`. The existing stage-dir vars (`BUILD` 224/180, `EM_DIR` 296/246,
  `DENSITY_DIR` 339/354, `EQUIL_DIR` 459, `PROD_DIR` 548/547, `HEAT_DIR` 289, `RELAX_DIR` 475)
  keep pointing at `${OUTDIR}/<stage>` — they now transparently resolve onto scratch for the
  symlinked ones. So **no mdrun/deffnm path edits** are needed.

### REMD engine — remove the now-redundant per-file symlink tricks
- Equil pre-run link (line ~509-510: `ln -sf "${SCRATCH_DIR}/equil_rep${idx}.xtc" "${REPDIR}/equil.xtc"`) → **delete**; `-deffnm equil` inside `equil/rep*/` already lands on scratch via the folder symlink.
- Prod pre-run link (line ~601-602: `ln -sf "${SCRATCH_DIR}/remd_rep${idx}.xtc" "${REPDIR}/remd.xtc"`) → **delete** (same reason).
- Density `-x "${SCRATCH_DIR}/..._density_seg${seg}.xtc"` (line ~410) → change to write under the
  (now symlinked) `${DENSITY_DIR}` via `-deffnm`, or just point `-x` at `${DENSITY_DIR}/...`.
- STEP 10 finalize symlink loops (lines 644-658) → **delete** (trajectories/ is gone).

### MD engine — same idea
- MD writes trajectories with `-x "${SCRATCH_DIR}/*.xtc"` directly (heat ~341, density ~414,
  relax ~523, prod ~599). Change these to write under the symlinked stage dirs
  (`${HEAT_DIR}/heat.xtc`, `${DENSITY_DIR}/..._seg.xtc`, `${RELAX_DIR}/relax.xtc`,
  `${PROD_DIR}/md.xtc`). STEP 8 finalize symlink loop (lines 612-618) → **delete**.

### New finalize step (both engines, last real step, after analysis)
- New `scripts/simulation/finalize_outputs.sh "$OUTDIR" "$SCRATCH_DIR"`:
  copy the **real** OUTDIR dirs/files (`analysis/ logs/ em/ build/ parameters.txt <base>_final*.pdb`)
  into `SCRATCH_DIR` at the same relative paths, so the scratch archive is self-contained.
  Bulk is already on scratch (written there during the run) → nothing to move. Echo the exact
  command first (per the print-before-run convention). Idempotent.
- Reached only on success (both engines run `set -euo pipefail`; put it after the analysis step).

### `parameters.txt` (REMD STEP 11 lines 676-704 / MD STEP 9)
- Already dual-written to OUTDIR + SCRATCH. Keep. Already records `Output directory` +
  `Scratch directory` (bidirectional pointer — satisfies "find the OUTDIR from a scratch dir").

### `run_analysis.sh`
- **Likely no change**: it reads `${OUTDIR}/prod/...` for TPR (real→resolves via symlink) and
  currently the trajectory from `${OUTDIR}/trajectories/<name>.xtc` (lines 64/69/76). Since
  `trajectories/` is dropped, update those 3 XTC lines to the stage path
  (`${OUTDIR}/prod/md.xtc`, `${OUTDIR}/prod/rep${REP}/remd.xtc`, `.../rest2.xtc`). Mode detection
  (`-f prod/…tpr`) and `-e "$XTC"` already follow symlinks. **REST2 engine uses the same pattern —
  update it in lockstep.**

### Cleanup / ERR trap
- The ERR trap `rm -rf "$SCRATCH_DIR"` on early failure (REMD 165-177 / MD 146-159) now deletes
  the *live run* on failure below `PRESERVE_FROM_STEP`. Since the bulk is on scratch from step 1,
  revisit the preserve thresholds so a failed run's scratch isn't nuked when it's the only copy of
  partial data worth inspecting. (Decision: likely lower/relax the threshold; confirm at impl time.)

## Risks / things to verify at implementation

- **Folder symlink + `-multidir` (REMD prod/equil):** confirm mdrun writes cleanly into
  `equil/rep*/` and `prod/rep*/` when those live under a symlinked parent (expected fine — it's a
  normal path). Smoke-test with the small example before a full run.
- **genion/solvate:** `build/` stays real on the pool → rename is same-fs. Do **not** symlink
  `build/` (would reintroduce the cross-fs rename gotcha).
- **`df` space check** (REMD 180 / MD) targets `$SCRATCH_DIR` — still correct (bulk goes there).
- **REST2 engine** (`REST2-gromacs.sbatch`) must get the same treatment to stay consistent.
- Update **CLAUDE.md**: document the folder-symlink contract, which stage dirs are symlinks vs
  real, why `build/` stays real (genion), and the laptop-rsync note (few folder-links to exclude).

## Test

Use the example run dir to exercise `finalize_outputs.sh` + the `run_analysis.sh` path change
without a full simulation, then a short REMD smoke run (small `example` PDB, few reps, ~50 ps) to
confirm the folder symlinks + `-multidir` + analysis-from-OUTDIR all work end-to-end.
