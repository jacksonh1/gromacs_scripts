# Plan: add the a99SB-disp force field

**Status:** OVERVIEW / not started (2026-10-08). Scoping only — what it would take to run
a99SB-disp in the T-REMD / MD / (eventually REST2) engines. No code written yet.

## What a99SB-disp is

a99SB-disp (Robustelli et al., *PNAS* 2018) is an Amber99SB\*-ILDN-Q protein force field
**co-parameterized with a modified 4-site TIP4P-D water** ("a99SBdisp water"), using c22
ions. It is designed to balance folded and intrinsically-disordered states. It is **not
bundled with GROMACS** — the GROMACS-format port is distributed from Paul Robustelli's
repository (<https://github.com/paulrobustelli/Force-Fields>). The protein side is
AMBER-family and uses **no CMAP**; the distinctive pieces are the 4-site water and modified
protein–water LJ (NBFIX-style `[ nonbond_params ]`).

## Work required (ordered by effort)

1. **Install the port (config only, mirrors CHARMM install Step 7).** Download the port,
   extract into `$GMXLIB` (`~/opt/gromacs/ff/`), keep the upstream directory name, optionally
   add an `FF_ALIASES` entry in `site_config.sh`. `WATER` is already a free job parameter
   (`validate_params.py` allowlist, no value restriction) — set it to whatever the port's
   `watermodels.dat` names the disp water.

2. **BLOCKER — 4-site water needs a real code change.** All three engines **hardcode**
   `gmx solvate -cs spc216.gro` (a 3-site box): `REMD-gromacs.sbatch:341`,
   `MD-gromacs.sbatch:290`, `REST2-gromacs.sbatch:377`. A 4-site water must be solvated from
   a 4-site box (`tip4p.gro` or the port's own water box). Fix = select the solvent box by
   water model via one shared helper, applied in all three engines. This is the main task.

3. **mdp nonbonded settings — verify, possibly a third branch.** The mdp gate is binary:
   `FF==charmm*` gets force-switch + `DispCorr=no`; everything else gets the AMBER branch
   (plain cutoff, `DispCorr=EnerPres`). a99SB-disp lands on the AMBER branch. Published
   settings fit mostly (cutoff **1.0 nm** → `CUTOFF_NM=1.0`, PME, fourierspacing 0.16 nm,
   LINCS h-bonds, 2 fs). **Open question: DispCorr.** TIP4P-D deliberately inflates C6
   dispersion; confirm from the port README whether Robustelli's protocol keeps analytic
   `DispCorr=EnerPres` or turns it off. If off, a third mdp case is needed (neither CHARMM
   nor plain-AMBER). See the "neither AMBER-like nor CHARMM-like" warning in
   `docs/FORCE_FIELDS.md`.

4. **Ions — verify, probably fine.** genion uses `-pname NA -nname CL`. The port recommends
   c22 ions; confirm those residue/atom names resolve in the port's `ions.itp` before relying
   on it.

5. **REST2 — run the scaling self-check before enabling.** No CMAP, so the STEP-1 `charmm*`
   guard does **not** auto-reject it. But the protein–water `[ nonbond_params ]` term means
   the `partial_tempering` λ=1.0-vs-0.5 diff check in `docs/FORCE_FIELDS.md` must be run to
   confirm every conformational solute term is scaled. Until that passes: T-REMD + MD only.

6. **Docs + test.** Add a row to the FORCE_FIELDS.md "Known status" table (record the
   measured nonbonded + REST2 result), note the port in installation README Step 7, then a
   throwaway `pdb2gmx` build followed by a short MD end-to-end (exercises the solvent-box fix
   and the mdp).

## Summary

One genuine code change (solvent-box selection for 4-site water, three engines), two things
to verify against the port README before any compute (DispCorr treatment, ion names), and the
REST2 scaling self-check before REST2 is allowed. Everything else is config + docs.

## References

- Port: <https://github.com/paulrobustelli/Force-Fields>
- Settings corroboration: a99SB-disp IDP study, <https://pmc.ncbi.nlm.nih.gov/articles/PMC10087257/>
- Related docs: `docs/FORCE_FIELDS.md` (FF resolution, mdp gate, REST2 check),
  `scripts/installation/README.md` Step 7 (installing a non-bundled FF).
