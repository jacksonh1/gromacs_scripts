#!/usr/bin/env python3
"""Validate per-job parameters for the GROMACS REMD / MD engines.

Two modes, selected by the presence of --check-keys:

  validate_params.py --engine {remd,md} --check-keys CONFIG
      Parse the assignment lines of a sourced shell config file and reject any
      key that is not a recognized parameter (catches typos like REPLICASS=48
      before the file is sourced and the typo silently uses a default).

  validate_params.py --engine {remd,md}
      Read the *resolved* parameters from the environment (after the engine's
      ${VAR:-default} resolution) and check that required ones are set, values
      are the right type, and ranges are sane. Covers params set either via a
      config file or via --export.

Pure standard library: this runs under the sbatch's system python3 at STEP 1,
before the analysis conda env is activated. On any problem it prints one or more
'[ERROR] ...' lines to stderr and exits 1 — never a traceback. On success it
prints '[OK] ...' and exits 0.
"""

import argparse
import os
import re
import sys

# ── Allowed job parameters per engine ────────────────────────────────────────
# The LHS names a user may set in a config file (or pass via --export). Keep in
# sync with the ${VAR:-default} reads in the matching *-gromacs.sbatch STEP 1.
_SHARED_JOB = {
    "PDB_IN", "OUTBASE", "OUTDIR", "FF", "WATER", "BOX_SHAPE", "BOX_BUFFER",
    "NEUTRALIZE", "SALT_MOLAR", "DT_PS", "CUTOFF_NM", "GAMMA_LN", "TOTAL_NS",
    "EQUI_NS", "NPT_SEG_STEPS", "NPT_MIN_SEG", "NPT_MAX_SEG", "NPT_TOL_REL",
    "PRESERVE_FROM_STEP", "SCRATCH_DIR", "SCRATCH_ROOT", "GMX",
}
_REMD_JOB = {"FORCE", "REPLICAS", "NTOMP_SERIAL", "T_MIN", "T_MAX", "TEMPS_LIST", "REPLEX_PS"}
_MD_JOB = {"T_SIM", "TRAJ_PS", "EQ_NPT_NS", "NTOMP"}

# site_config.sh settings a job config may legitimately override (so they are not
# flagged as typos), even though they are normally set once per cluster.
_SITE_OVERRIDES = {
    "GMXRC", "PLUMED_SH", "CUDA_MODULE", "OPENMPI_MODULE", "CONDA_MODULE",
    "GROMD_ENV", "GROMACS_SCRIPTS_DIR",
}

ALLOWED = {
    "remd": _SHARED_JOB | _REMD_JOB | _SITE_OVERRIDES,
    "md": _SHARED_JOB | _MD_JOB | _SITE_OVERRIDES,
}

_ASSIGN = re.compile(r"^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=")


def check_keys(path, engine):
    """Return a list of [ERROR] strings for unrecognized assignment keys."""
    allowed = ALLOWED[engine]
    errs = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            stripped = line.lstrip()
            if not stripped or stripped.startswith("#"):
                continue
            m = _ASSIGN.match(line)
            if not m:
                continue
            key = m.group(1)
            if key not in allowed:
                errs.append(
                    f"[ERROR] Unknown config key '{key}' on line {lineno} of {path} "
                    f"(typo? not a recognized {engine.upper()} parameter)"
                )
    return errs


def validate_values(engine):
    """Return a list of [ERROR] strings for missing/wrong-type/out-of-range values."""
    env = os.environ
    errs = []

    def as_num(name, raw, integer):
        try:
            return int(raw) if integer else float(raw)
        except ValueError:
            errs.append(f"[ERROR] {name}='{raw}' is not {'an integer' if integer else 'a number'}")
            return None

    def positive(name, *, integer=False, allow_zero=False):
        """Validate an optional numeric env var is > 0 (or >= 0). Returns the value or None."""
        raw = env.get(name, "")
        if raw == "":
            return None  # unset → the engine default applies
        v = as_num(name, raw, integer)
        if v is None:
            return None
        if v < 0 or (v == 0 and not allow_zero):
            errs.append(f"[ERROR] {name} must be {'>= 0' if allow_zero else '> 0'} (got '{raw}')")
        return v

    def flag01(name):
        raw = env.get(name, "")
        if raw not in ("", "0", "1"):
            errs.append(f"[ERROR] {name} must be 0 or 1 (got '{raw}')")

    # ── required ──
    pdb = env.get("PDB_IN", "")
    if pdb == "":
        errs.append("[ERROR] PDB_IN is required (no default) — set the input structure")
    elif not os.path.isfile(pdb):
        errs.append(f"[ERROR] PDB_IN file not found: {pdb}")

    # ── shared numeric ──
    for name in ("DT_PS", "CUTOFF_NM", "GAMMA_LN", "BOX_BUFFER", "TOTAL_NS", "NPT_TOL_REL"):
        positive(name)
    positive("EQUI_NS", allow_zero=True)
    nmin = positive("NPT_MIN_SEG", integer=True)
    nmax = positive("NPT_MAX_SEG", integer=True)
    positive("NPT_SEG_STEPS", integer=True)
    if nmin is not None and nmax is not None and nmin > nmax:
        errs.append(f"[ERROR] NPT_MIN_SEG ({nmin}) must be <= NPT_MAX_SEG ({nmax})")
    salt = env.get("SALT_MOLAR", "")
    if salt != "":
        v = as_num("SALT_MOLAR", salt, integer=False)
        if v is not None and v < 0:
            errs.append(f"[ERROR] SALT_MOLAR must be >= 0 (got '{salt}')")
    flag01("NEUTRALIZE")

    # ── engine-specific ──
    if engine == "remd":
        flag01("FORCE")
        positive("REPLEX_PS")
        tmin = positive("T_MIN")
        tmax = positive("T_MAX")
        if tmin is not None and tmax is not None and tmax <= tmin:
            errs.append(f"[ERROR] T_MAX ({tmax}) must be > T_MIN ({tmin})")
        rep = positive("REPLICAS", integer=True)
        if rep is not None and rep < 2:
            errs.append(f"[ERROR] REPLICAS must be >= 2 (got {rep})")
        # A config-supplied REPLICAS is silently overridden by SLURM's -n (ntasks).
        cfg_rep = env.get("CFG_REPLICAS", "")
        ntasks = env.get("SLURM_NTASKS", "")
        if cfg_rep and ntasks and cfg_rep != ntasks:
            errs.append(
                f"[ERROR] config REPLICAS ({cfg_rep}) != SLURM ntasks ({ntasks}); the run "
                f"uses ntasks, so the config value is ignored — submit with -n {cfg_rep} "
                f"or set REPLICAS={ntasks}"
            )
    else:  # md
        positive("T_SIM")
        positive("TRAJ_PS")
        positive("EQ_NPT_NS", allow_zero=True)
        nt = positive("NTOMP", integer=True)
        if nt is not None and nt < 1:
            errs.append(f"[ERROR] NTOMP must be >= 1 (got {nt})")

    return errs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", required=True, choices=["remd", "md"])
    ap.add_argument("--check-keys", metavar="CONFIG",
                    help="validate the assignment keys in this config file instead of env values")
    args = ap.parse_args()

    if args.check_keys is not None:
        if not os.path.isfile(args.check_keys):
            print(f"[ERROR] config file not found: {args.check_keys}", file=sys.stderr)
            sys.exit(1)
        errs = check_keys(args.check_keys, args.engine)
        label = f"config keys in {args.check_keys}"
    else:
        errs = validate_values(args.engine)
        label = f"{args.engine.upper()} parameters"

    if errs:
        for e in errs:
            print(e, file=sys.stderr)
        sys.exit(1)
    print(f"[OK] {label} validated")


if __name__ == "__main__":
    main()
