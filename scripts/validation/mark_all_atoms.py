#!/usr/bin/env python3
"""VALIDATION ONLY: mark EVERY atom of a processed topology as hot (append `_` to its type).

With everything hot, `partial_tempering <lambda>` must scale every non-bonded and torsion
energy term by exactly lambda while leaving bonds and angles alone — a clean ratio test of
the scaling (validate_rest2_scaling.sh, check (b)). Never use this for a production run:
REST2 heats the solute only (mark_hot_region.py).

Usage: mark_all_atoms.py processed.top > allhot.top
"""
import re
import sys

section = ""
marked = 0
with open(sys.argv[1]) as fh:
    for line in fh:
        text = line.split(";", 1)[0].strip()
        m = re.fullmatch(r"\[\s*(\w+)\s*\]", text)
        if m:
            section = m.group(1).lower()
        elif section == "atoms" and text:
            f = line.split()
            if not f[1].endswith("_"):
                f[1] += "_"
            marked += 1
            line = " ".join(f) + "\n"
        sys.stdout.write(line)
if marked == 0:
    raise SystemExit("[ERROR] mark_all_atoms: no [ atoms ] lines found")
print(f"[OK] mark_all_atoms: marked {marked} atoms", file=sys.stderr)
