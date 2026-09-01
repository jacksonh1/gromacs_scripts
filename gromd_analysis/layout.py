# ─────────────────────────────────────────────────────────────────────────────
# layout.py — parse a finished job directory into a typed JobDir
# ─────────────────────────────────────────────────────────────────────────────
# One place that knows what a finished MD / T-REMD / REST2 job directory looks
# like. Everything downstream — run_analysis.sh, notebooks, one-off analyses —
# takes a JobDir instead of re-deriving paths from a bare OUTDIR.
#
# This is a parser, not a validator: JobDir.detect() either returns a JobDir
# whose paths are known to exist, or raises. Holding a JobDir *is* the proof
# that the directory was a real job, and it carries which kind it was.
#
# Usage (library):
#   from gromd_analysis.layout import JobDir
#   job = JobDir.detect(Path(outdir), rep="000")
#
# Usage (shell, via the gromd-layout entry point):
#   eval "$(gromd-layout OUTDIR [REP])"
#   # sets MODE TPR XTC PREFIX NCHAINS ANALYSIS_DIR EM_GRO EM_TPR BUILD_DIR
#
# Informational messages go to stderr — stdout is eval'd by the caller.
# ─────────────────────────────────────────────────────────────────────────────

import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from gromd_analysis.chains import parse_molecules

Mode = Literal["MD", "REMD", "REST2"]

# Production basenames are method-named, not ensemble-named; the analysis layer
# keys off them to tell the engines apart. (stage, basename) per mode:
#   MD    prod/md.tpr          REMD  prod/rep<REP>/remd.tpr
#   REST2 prod/rep<REP>/rest2.tpr  — rep000 is lambda=1, the physical T_MIN ensemble,
#         so its analysis is identical to REMD's slot analysis; only the stem differs.
_MD_BASENAME = "md"
_PER_REPLICA_BASENAMES: dict[str, Mode] = {"remd": "REMD", "rest2": "REST2"}


class JobLayoutError(Exception):
    """The directory is not a finished job, or its run input/trajectory is missing."""


@dataclass(frozen=True)
class JobDir:
    """A finished job directory, resolved. Every path here exists.

    prefix is a path *stem*, not a file: analysis outputs are written as
    <prefix>_rmsd.xvg, <prefix>_stripped_aligned.xtc, and so on.
    """

    root: Path
    mode: Mode
    tpr: Path
    xtc: Path
    prefix: Path
    n_chains: int
    analysis_dir: Path
    build_dir: Path
    em_gro: Path
    em_tpr: Path

    @property
    def multichain(self) -> bool:
        return self.n_chains > 1

    @classmethod
    def detect(cls, outdir: Path, rep: str = "000") -> "JobDir":
        root = Path(outdir).resolve()
        if not root.is_dir():
            raise JobLayoutError(f"Not a directory: {outdir}")

        analysis_dir = root / "analysis"
        prod = root / "prod"

        if (prod / f"{_MD_BASENAME}.tpr").is_file():
            mode: Mode = "MD"
            stage = prod
            basename = _MD_BASENAME
            prefix = analysis_dir / _MD_BASENAME
        else:
            stage = prod / f"rep{rep}"
            for basename, candidate_mode in _PER_REPLICA_BASENAMES.items():
                if (stage / f"{basename}.tpr").is_file():
                    mode = candidate_mode
                    prefix = analysis_dir / f"{basename}_rep{rep}"
                    break
            else:
                raise JobLayoutError(
                    f"Cannot find prod/md.tpr, prod/rep{rep}/remd.tpr, or "
                    f"prod/rep{rep}/rest2.tpr under {root}\n"
                    f"        Is this a finished MD, T-REMD, or REST2 job directory?"
                )

        tpr = stage / f"{basename}.tpr"
        xtc = stage / f"{basename}.xtc"

        # Jobs run before the output restructure collected trajectories in a top-level
        # trajectories/ dir. That dir is gone from the current pipeline, but old job
        # dirs (including the shipped example) still have it — accept it explicitly and
        # announce it, rather than treating it as a search path.
        if not xtc.exists() and (root / "trajectories" / xtc.name).exists():
            xtc = root / "trajectories" / xtc.name
            print(f"[INFO] Legacy layout: trajectory taken from trajectories/ → {xtc}",
                  file=sys.stderr)

        if not tpr.is_file():
            raise JobLayoutError(f"Run input not found: {tpr}")
        # .exists() follows symlinks: under SYMLINK_BULK=1 prod/ is a symlink into
        # scratch, so the .xtc resolves through it. A missing file means scratch was purged.
        if not xtc.exists():
            raise JobLayoutError(f"Trajectory not found (scratch purged?): {xtc}")

        build_dir = root / "build"
        return cls(
            root=root,
            mode=mode,
            tpr=tpr,
            xtc=xtc,
            prefix=prefix,
            n_chains=_count_protein_chains(build_dir),
            analysis_dir=analysis_dir,
            build_dir=build_dir,
            em_gro=root / "em" / "em.gro",
            em_tpr=root / "em" / "em.tpr",
        )


def _count_protein_chains(build_dir: Path) -> int:
    """Protein molecule count from the system topology's [ molecules ] section.

    >1 means the complex must be kept in one periodic image (the multichain_*
    pipeline). Falls back to 1 when the topology is unreadable — the same
    assumption the pipeline has always made, but said out loud, because picking
    the single-chain path for a real complex silently inflates RMSD.
    """
    tops = sorted(build_dir.glob("*.top"))
    if not tops:
        print(f"[WARN] No .top in {build_dir}; assuming a single protein chain.",
              file=sys.stderr)
        return 1

    total = sum(count for moltype, count in parse_molecules(tops[0])
                if moltype.startswith("Protein"))
    if total == 0:
        print(f"[WARN] No Protein rows in {tops[0]}; assuming a single protein chain.",
              file=sys.stderr)
        return 1
    return total


def main() -> None:
    """Emit a JobDir as shell assignments for `eval`. All chatter goes to stderr."""
    if len(sys.argv) not in (2, 3):
        sys.exit("Usage: gromd-layout OUTDIR [REP]")
    outdir = Path(sys.argv[1])
    rep = sys.argv[2] if len(sys.argv) == 3 else "000"

    try:
        job = JobDir.detect(outdir, rep)
    except JobLayoutError as exc:
        sys.exit(f"[ERROR] {exc}")

    for var, value in [
        ("MODE", job.mode),
        ("TPR", job.tpr),
        ("XTC", job.xtc),
        ("PREFIX", job.prefix),
        ("NCHAINS", job.n_chains),
        ("ANALYSIS_DIR", job.analysis_dir),
        ("EM_GRO", job.em_gro),
        ("EM_TPR", job.em_tpr),
        ("BUILD_DIR", job.build_dir),
    ]:
        print(f"{var}={shlex.quote(str(value))}")


if __name__ == "__main__":
    main()
