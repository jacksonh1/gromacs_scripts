#!/usr/bin/env python3
"""Decide whether the iterative NPT density equilibration has reached a plateau.

Usage:
    density_converged.py --tol TOL --min-seg N V1 [V2 ...]

V1..Vn are the per-segment average volumes (nm^3), oldest first. Prints ONE line
to stdout: the verdict ('1' converged / '0' not), a space, then a human-readable
reason. The caller splits it -- verdict into its `converged` test, reason into the
job log next to the per-segment volumes. Exit status is 0 either way; a non-zero
exit means the arguments were bad, not that the density has not converged.

Both halves go to stdout on purpose. An earlier version sent the reason to stderr,
which split the two halves of one decision across the SLURM .out and .err files.

WHY A SLOPE TEST, NOT A CONSECUTIVE-SEGMENT DIFFERENCE
------------------------------------------------------
The engines used to stop as soon as |V_n - V_{n-1}| / V_{n-1} <= tol. For a
solvated box the segment-to-segment noise is well below that threshold
(sigma_V/V ~ sqrt(kT*kappa/V) ~ 2e-3, and the standard error of a segment mean
is smaller still), so the old test passes readily -- including while the box is
still shrinking steadily. A sustained drift of 0.4 % per segment clears a 0.5 %
consecutive-difference threshold on every single comparison while the box
contracts several percent overall.

This fits a least-squares line through the trailing window instead and asks how
much systematic drift that slope accounts for ACROSS the whole window. Random
scatter has ~zero slope and passes; a steady drift does not, however small each
individual step is. Scatter itself is deliberately NOT part of the test -- it is
the physical volume fluctuation of an NPT box, not a sign of non-convergence.

Pure standard library: this runs under the sbatch's system python3, before the
analysis conda env is activated.
"""

import argparse
import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class PlateauTest:
    """Result of the trailing-window drift test."""

    converged: bool
    window: int          # segments actually used
    drift_rel: float     # |slope| * (window - 1) / mean, i.e. fractional drift across the window
    mean_volume: float   # nm^3, over the window


def drift_over_window(volumes: list[float], window: int) -> PlateauTest:
    """Fractional volume drift the best-fit line accounts for across the last `window` segments.

    ASSUMES: volumes are per-segment averages in chronological order and the
    segments are of equal duration, so the segment index is a valid time axis.
    """
    assert window >= 2, "a slope needs at least two points"
    assert len(volumes) >= window, "caller must supply at least `window` volumes"

    tail = volumes[-window:]
    mean_v = sum(tail) / window
    assert mean_v > 0, "volumes must be positive"

    # OLS slope of V against segment index (x = 0..window-1, so x is centred by hand).
    mean_x = (window - 1) / 2.0
    sxy = sum((x - mean_x) * (v - mean_v) for x, v in enumerate(tail))
    sxx = sum((x - mean_x) ** 2 for x in range(window))
    slope = sxy / sxx

    drift_rel = abs(slope) * (window - 1) / mean_v
    return PlateauTest(converged=False, window=window, drift_rel=drift_rel, mean_volume=mean_v)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol", type=float, required=True,
                    help="max fractional drift across the window (e.g. 0.005)")
    ap.add_argument("--min-seg", type=int, required=True,
                    help="trailing window length in segments; also the earliest segment "
                         "at which convergence may be declared")
    ap.add_argument("volumes", nargs="+", type=float,
                    help="per-segment average volumes (nm^3), oldest first")
    args = ap.parse_args()

    if args.tol <= 0:
        sys.exit(f"[ERROR] --tol must be > 0 (got {args.tol})")
    if args.min_seg < 2:
        sys.exit(f"[ERROR] --min-seg must be >= 2 (got {args.min_seg})")

    if len(args.volumes) < args.min_seg:
        print(f"0 {len(args.volumes)}/{args.min_seg} segments — too few to test")
        return

    test = drift_over_window(args.volumes, args.min_seg)
    converged = test.drift_rel <= args.tol
    print(f"{1 if converged else 0} "
          f"drift over last {test.window} segments = "
          f"{test.drift_rel:.2%} of {test.mean_volume:.2f} nm^3 "
          f"(tol {args.tol:.2%}) — {'converged' if converged else 'still drifting'}")


if __name__ == "__main__":
    main()
