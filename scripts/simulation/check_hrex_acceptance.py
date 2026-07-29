#!/usr/bin/env python3
"""Fail-loud gate: does a GROMACS replica-exchange log show ANY accepted exchanges?

Parses the 'Replica exchange statistics' block of a replica log and looks at the
'average probabilities' row. Exits 0 if at least one pair has probability > 0,
exits 1 (with an [ERROR]) if every pair is 0.00 — the signature of a broken hrex
build (e.g. the GROMACS-2024 PLUMED patch), where REST2 runs but silently never
exchanges. The REST2 engine calls this on a short pre-flight run so a misconfigured
build can never masquerade as a working REST2 production run.

stdlib only (runs under the sbatch's system python3).

Usage: check_hrex_acceptance.py <replica.log>
"""
import re
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: check_hrex_acceptance.py <replica.log>")
    log = sys.argv[1]
    try:
        lines = open(log).read().splitlines()
    except OSError as e:
        print(f"[ERROR] cannot read log: {e}", file=sys.stderr)
        sys.exit(1)

    # Find the LAST 'average probabilities:' block (end-of-run cumulative stats).
    idx = [i for i, l in enumerate(lines) if "average probabilities" in l]
    if not idx:
        print(f"[ERROR] no 'average probabilities' block in {log} — did replica exchange run?",
              file=sys.stderr)
        sys.exit(1)
    start = idx[-1]

    # Format:
    #   Repl  average probabilities:
    #   Repl     0    1    2 ...        <- replica indices
    #   Repl      .74  .65 ...          <- the probabilities
    # The values row is the first 'Repl' line after the header whose tokens parse as floats.
    probs = None
    for l in lines[start + 1:start + 5]:
        toks = l.split()
        if not toks or toks[0] != "Repl":
            break
        vals = toks[1:]
        # index row is integers (0,1,2,...); values row has decimals like .74 / 0.65
        if vals and all(re.fullmatch(r"\d+", v) for v in vals):
            continue  # index row, skip
        try:
            probs = [float(v) for v in vals]
            break
        except ValueError:
            continue

    if probs is None:
        print(f"[ERROR] could not parse the average-probabilities row in {log}", file=sys.stderr)
        sys.exit(1)

    nonzero = [p for p in probs if p > 0.0]
    if not nonzero:
        print(f"[ERROR] hrex acceptance is ZERO for all {len(probs)} pairs in {log}.\n"
              f"        REST2 exchange is not working — almost certainly the wrong GROMACS build\n"
              f"        (hrex is broken on the 2024.3-plumed build). Check REST2_GMXRC.",
              file=sys.stderr)
        sys.exit(1)

    print(f"[OK] hrex acceptance nonzero: {len(nonzero)}/{len(probs)} pairs exchanging "
          f"(avg prob {sum(probs)/len(probs):.2f})")
    sys.exit(0)


if __name__ == "__main__":
    main()
