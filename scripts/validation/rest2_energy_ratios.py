#!/usr/bin/env python3
"""VALIDATION ONLY: compare per-term rerun energies of a scaled REST2 topology to the unscaled one.

Usage:
    rest2_energy_ratios.py REF.xvg TEST.xvg EXPECT [--tol 1e-4]

REF/TEST are `gmx energy` outputs over the SAME frames with the SAME term selection (the
legend names must match). EXPECT is comma-separated `term=ratio` pairs: `term` is the exact
`gmx energy` legend name and `ratio` a number or `sqrt(x)`, e.g.
"Bond=1,Angle=1,Proper Dih.=0.5,LJ-SR:Protein-non-Protein=sqrt(0.5)", or `*` to report a
term without checking it (a total that mixes hot and cold pairs has no single ratio; its
energy-group parts carry the check). Every legend term
must have an EXPECT entry and every EXPECT entry must match a term, so a term that
silently vanished (or appeared) fails instead of going unchecked. The ratio is checked
frame by frame: TEST = ratio * REF within tol (relative) on every frame.
Exit 1 on any mismatch.
"""
import argparse
import math
import re
import sys


def read_xvg(path):
    names, rows = [], []
    for ln in open(path):
        m = re.match(r'@ s\d+ legend "(.*)"', ln)
        if m:
            names.append(m.group(1))
        elif ln[0] not in "#@":
            rows.append([float(v) for v in ln.split()[1:]])
    if not rows or len(rows[0]) != len(names):
        raise SystemExit(f"[ERROR] {path}: {len(names)} legends vs {len(rows[0]) if rows else 0} columns")
    return names, rows


def parse_expect(spec):
    out = []
    for item in spec.split(","):
        term, _, val = item.partition("=")
        val = val.strip()
        m = re.fullmatch(r"sqrt\((.*)\)", val)
        out.append((term.strip(), None if val == "*" else math.sqrt(float(m.group(1))) if m else float(val)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("test")
    ap.add_argument("expect")
    ap.add_argument("--tol", type=float, default=1e-4)
    a = ap.parse_args()
    rn, rr = read_xvg(a.ref)
    tn, tr = read_xvg(a.test)
    if rn != tn or len(rr) != len(tr):
        raise SystemExit(f"[ERROR] REF and TEST differ in terms or frames:\n  {rn}\n  {tn}")
    expect = parse_expect(a.expect)
    used = set()
    bad = 0
    print(f"{'term':34s} {'expected':>9s} {'min ratio':>10s} {'max ratio':>10s}  verdict")
    for j, name in enumerate(rn):
        hits = [(t, v) for t, v in expect if t == name]
        if len(hits) != 1:
            print(f"{name:34s} has {len(hits)} EXPECT entries (need exactly 1)")
            bad += 1
            continue
        term, want = hits[0]
        used.add(term)
        ratios, ok = [], True
        if want is None:
            rs = [t[j] / r[j] for r, t in zip(rr, tr, strict=True) if abs(r[j]) > 1e-9]
            print(f"{name:34s} {'*':>9s} {min(rs):10.6f} {max(rs):10.6f}  (not checked)")
            continue
        for ref_row, test_row in zip(rr, tr, strict=True):
            r, t = ref_row[j], test_row[j]
            if abs(r) < 1e-9:
                ok &= abs(t) < 1e-6
                continue
            ratios.append(t / r)
            ok &= abs(t - want * r) <= a.tol * abs(r) + 1e-6
        lo, hi = (min(ratios), max(ratios)) if ratios else (float("nan"),) * 2
        print(f"{name:34s} {want:9.6f} {lo:10.6f} {hi:10.6f}  {'ok' if ok else 'FAIL'}")
        bad += not ok
    for term, _ in expect:
        if term not in used:
            print(f"EXPECT entry '{term}' matched no energy term")
            bad += 1
    if bad:
        print(f"[FAIL] {bad} term(s) off", file=sys.stderr)
        return 1
    print(f"[OK] all {len(rn)} terms at the expected ratio over {len(rr)} frame(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
