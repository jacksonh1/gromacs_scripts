#!/usr/bin/env python3
"""Verify that `plumed partial_tempering` scaled a REST2 topology the way REST2 requires.

Usage:
    check_rest2_topology.py COLD.top HOT.top --lambda-cold 1.0 --lambda-hot 0.5

COLD/HOT are two outputs of `partial_tempering` made from the SAME marked topology
(STEP 7: topol_rep000.top and the replica with the lowest lambda). Hot atoms are the ones
whose atom type carries the `_` suffix (mark_hot_region.py). Checked, per the REST2
definition and the partial_tempering source:

  [ atoms ]           hot charge = cold charge * sqrt(r)            (r = lambda_hot/lambda_cold)
  [ dihedrals ]       force constants * sqrt(r) per hot atom in positions 1 and 4
  [ atomtypes ]       TYPE_ epsilon = lambda * epsilon(TYPE), sigma unchanged   (each file)
  [ nonbond_params ]  / [ pairtypes ]: epsilon * lambda^(n_hot/2), sigma unchanged (each file)
  bonds, angles, and every other line: identical between COLD and HOT

and it fails on anything partial_tempering cannot scale: a section such as [ cmap ] /
[ cmaptypes ], or a [ dihedrals ] function type outside {1,2,3,4,5,9}. Those would leave
part of the solute Hamiltonian at full strength while every other check still passes
(knowledgebase/GOTCHAS.md, "partial_tempering does NOT scale CHARMM CMAP").

Exit 0 and print a summary on success; exit 1 listing the offending lines otherwise.
Stdlib-only, like the other STEP helpers.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass

SUFFIX = "_"

# Sections partial_tempering never rewrites but that carry solute conformational energy.
RED_FLAG_SECTIONS = frozenset({
    "cmap", "cmaptypes", "pairs_nb",
    "polarization", "thole_polarization", "water_polarization",
    "dihedral_restraints", "distance_restraints", "orientation_restraints",
    "angle_restraints", "angle_restraints_z",
})
# [ dihedrals ] function type -> 0-based field indices partial_tempering scales
# (fields 0-3 atoms, 4 funct; see partial_tempering.sh "DIHEDRALS").
DIHEDRAL_SCALED = {1: (6,), 4: (6,), 9: (6,), 2: (6,), 3: (5, 6, 7, 8, 9, 10), 5: (5, 6, 7, 8)}
DIHEDRAL_NFIELDS = {1: 8, 4: 8, 9: 8, 2: 7, 3: 11, 5: 9}

REL_TOL = 2e-5   # partial_tempering reprints numbers with awk's ~6 significant digits
ABS_TOL = 2e-6


@dataclass(frozen=True)
class Line:
    lineno: int
    section: str
    fields: tuple[str, ...]


def parse(path: str) -> list[Line]:
    out: list[Line] = []
    section = ""
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            text = raw.split(";", 1)[0].strip()
            if not text:
                continue
            if text.startswith("#"):
                raise SystemExit(f"[ERROR] {path}:{lineno}: preprocessor line in a processed topology: {text}")
            m = re.fullmatch(r"\[\s*(\w+)\s*\]", text)
            if m:
                section = m.group(1).lower()
                continue
            out.append(Line(lineno, section, tuple(text.split())))
    return out


def close(a: float, b: float) -> bool:
    return abs(a - b) <= ABS_TOL + REL_TOL * max(abs(a), abs(b))


class Report:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.counts: dict[str, int] = {}

    def fail(self, msg: str) -> None:
        self.errors.append(msg)

    def count(self, key: str) -> None:
        self.counts[key] = self.counts.get(key, 0) + 1


def check_self_scaling(lines: list[Line], lam: float, tag: str, rep: Report) -> None:
    """Within one file: TYPE_ atomtypes and suffixed nonbond_params/pairtypes vs their originals."""
    combrule = next(int(ln.fields[1]) for ln in lines if ln.section == "defaults")
    types = {ln.fields[0]: ln for ln in lines if ln.section == "atomtypes"}
    for name, ln in types.items():
        if not name.endswith(SUFFIX):
            continue
        orig = types.get(name[: -len(SUFFIX)])
        if orig is None:
            rep.fail(f"{tag}:{ln.lineno}: hot atomtype {name} has no unscaled original")
            continue
        s_h, e_h = float(ln.fields[-2]), float(ln.fields[-1])
        s_o, e_o = float(orig.fields[-2]), float(orig.fields[-1])
        want_s = s_o * lam if combrule == 1 else s_o
        if not (close(e_h, e_o * lam) and close(s_h, want_s)):
            rep.fail(f"{tag}:{ln.lineno}: atomtype {name} sigma/eps {s_h} {e_h} != expected {want_s} {e_o * lam}")
        rep.count(f"{tag} atomtypes scaled")
    for sec in ("nonbond_params", "pairtypes"):
        table = {(ln.fields[0], ln.fields[1]): ln for ln in lines if ln.section == sec}
        for (a, b), ln in table.items():
            n_hot = a.endswith(SUFFIX) + b.endswith(SUFFIX)
            if n_hot == 0:
                continue
            orig = table.get((a.removesuffix(SUFFIX), b.removesuffix(SUFFIX)))
            if orig is None:
                rep.fail(f"{tag}:{ln.lineno}: [ {sec} ] {a} {b} has no unscaled original")
                continue
            f = lam ** (n_hot / 2)
            c6_h, c12_h = float(ln.fields[3]), float(ln.fields[4])
            c6_o, c12_o = float(orig.fields[3]), float(orig.fields[4])
            want_c6 = c6_o * f if combrule == 1 else c6_o
            if not (close(c12_h, c12_o * f) and close(c6_h, want_c6)):
                rep.fail(f"{tag}:{ln.lineno}: [ {sec} ] {a} {b} {c6_h} {c12_h} != expected {want_c6} {c12_o * f}")
            rep.count(f"{tag} {sec} scaled")


def check_pair(cold: list[Line], hot: list[Line], r: float, rep: Report) -> None:
    """Line-by-line COLD vs HOT: only hot charges and dihedral force constants may differ."""
    if len(cold) != len(hot):
        rep.fail(f"COLD has {len(cold)} data lines, HOT has {len(hot)}: not made from the same marked topology")
        return
    hot_atoms: set[str] = set()
    for c, h in zip(cold, hot, strict=True):
        if c.section != h.section or len(c.fields) != len(h.fields):
            rep.fail(f"structure differs at COLD:{c.lineno} / HOT:{h.lineno}")
            return
        sec = c.section
        if sec == "moleculetype":
            hot_atoms = set()
        if sec in RED_FLAG_SECTIONS:
            if not rep.counts.get(f"red flag {sec}"):
                rep.fail(f"HOT:{h.lineno}: [ {sec} ] is present and partial_tempering does not scale it")
            rep.count(f"red flag {sec}")
            continue
        if sec in ("atomtypes", "nonbond_params", "pairtypes"):
            continue  # scaled relative to the same file's originals: check_self_scaling
        if sec == "atoms":
            is_hot = h.fields[1].endswith(SUFFIX)
            if is_hot:
                hot_atoms.add(h.fields[0])
                q_c, q_h = float(c.fields[6]), float(h.fields[6])
                if not close(q_h, q_c * math.sqrt(r)):
                    rep.fail(f"HOT:{h.lineno}: hot atom {h.fields[0]} charge {q_h} != {q_c} * sqrt({r})")
                if c.fields[:6] != h.fields[:6] or c.fields[7:] != h.fields[7:]:
                    rep.fail(f"HOT:{h.lineno}: hot atom {h.fields[0]} changed outside the charge column")
                rep.count("hot atom charges scaled")
            elif c.fields != h.fields:
                rep.fail(f"HOT:{h.lineno}: cold atom {h.fields[0]} changed")
            continue
        if sec == "dihedrals":
            funct = int(h.fields[4])
            if funct not in DIHEDRAL_SCALED:
                rep.fail(f"HOT:{h.lineno}: [ dihedrals ] funct {funct} is not scaled by partial_tempering")
                continue
            if len(h.fields) != DIHEDRAL_NFIELDS[funct]:
                rep.fail(f"HOT:{h.lineno}: [ dihedrals ] funct {funct} without inline parameters")
                continue
            n_hot = (h.fields[0] in hot_atoms) + (h.fields[3] in hot_atoms)
            f = r ** (n_hot / 2)
            scaled = DIHEDRAL_SCALED[funct]
            for i, (a, b) in enumerate(zip(c.fields, h.fields, strict=True)):
                if i in scaled:
                    if not close(float(b), float(a) * f):
                        rep.fail(f"HOT:{h.lineno}: dihedral parameter {b} != {a} * {f:.6g}")
                elif a != b and not close(float(a), float(b)):
                    rep.fail(f"HOT:{h.lineno}: dihedral field {i} changed ({a} -> {b})")
            rep.count(f"dihedrals with {n_hot} hot end atom(s)")
            continue
        if c.fields != h.fields:
            rep.fail(f"HOT:{h.lineno}: [ {sec} ] line changed but REST2 must not scale it: {' '.join(h.fields)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cold")
    ap.add_argument("hot")
    ap.add_argument("--lambda-cold", type=float, required=True)
    ap.add_argument("--lambda-hot", type=float, required=True)
    a = ap.parse_args()
    # lambda-hot == lambda-cold is legitimate (the identical-lambda exchange sanity run):
    # the checks then demand hot == cold, which is what partial_tempering must produce.
    if not (0 < a.lambda_hot <= a.lambda_cold <= 1):
        raise SystemExit(f"[ERROR] need 0 < lambda-hot <= lambda-cold <= 1, got {a.lambda_hot} {a.lambda_cold}")
    cold, hot = parse(a.cold), parse(a.hot)
    rep = Report()
    check_self_scaling(cold, a.lambda_cold, "COLD", rep)
    check_self_scaling(hot, a.lambda_hot, "HOT", rep)
    check_pair(cold, hot, a.lambda_hot / a.lambda_cold, rep)
    required = ("hot atom charges scaled", "HOT atomtypes scaled")
    for key in required:
        if not rep.counts.get(key):
            rep.fail(f"nothing counted for '{key}' - is the hot region marked?")
    if not any(k.startswith("dihedrals with") and k != "dihedrals with 0 hot end atom(s)" for k in rep.counts):
        rep.fail("no dihedral was scaled - is the hot region marked?")
    if rep.errors:
        print(f"[ERROR] check_rest2_topology: {len(rep.errors)} problem(s) in {a.hot} vs {a.cold}:", file=sys.stderr)
        for e in rep.errors[:20]:
            print(f"  {e}", file=sys.stderr)
        if len(rep.errors) > 20:
            print(f"  ... and {len(rep.errors) - 20} more", file=sys.stderr)
        return 1
    summary = ", ".join(f"{k}: {v}" for k, v in sorted(rep.counts.items()))
    print(f"[OK] check_rest2_topology (lambda {a.lambda_cold} -> {a.lambda_hot}): {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
