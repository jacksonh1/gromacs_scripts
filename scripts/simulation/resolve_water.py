#!/usr/bin/env python3
"""Resolve a WATER name against a force field: the pdb2gmx arguments and the solvent box.

Usage:
    resolve_water.py FF WATER DIR [DIR ...]

DIRs are searched in order for FF.ff (pass the GMXLIB entries, then $GMXDATA/top — the
order pdb2gmx itself uses). Prints ONE line for the engine to split:

    <pdb2gmx -water arg> <stdin answer or -> <n sites> <solvent box>

Why this exists: `pdb2gmx -water` accepts only GROMACS's built-in names (spc, spce, tip3p,
tip4p, tip4pew, tip5p, tips3p), so a force-field-specific water such as a99SB-disp's
`a99SBdisp_water` must be chosen through `-water select` with its menu index from the FF's
watermodels.dat. And `solvate -cs` must use a box with the model's site count — a 3-site
spc216.gro box under a 4-site topology is a fatal grompp mismatch. A box shipped in the FF
dir as WATER.gro wins (its atom names match the FF's .itp); otherwise the site count picks
GROMACS's spc216.gro / tip4p.gro / tip5p.gro. Stdlib-only, like the other STEP helpers.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

PDB2GMX_BUILTIN = frozenset({"spc", "spce", "tip3p", "tip4p", "tip4pew", "tip5p", "tips3p"})
BOX_BY_SITES = {3: "spc216.gro", 4: "tip4p.gro", 5: "tip5p.gro"}


def find_ffdir(ff: str, dirs: list[str]) -> Path:
    for d in dirs:
        cand = Path(d) / f"{ff}.ff"
        if cand.is_dir():
            return cand
    raise SystemExit(f"[ERROR] resolve_water: {ff}.ff not found in: {' '.join(dirs)}")


def water_menu(ffdir: Path) -> list[str]:
    """Model names in watermodels.dat order: pdb2gmx's `-water select` menu is 1-based on this."""
    path = ffdir / "watermodels.dat"
    if not path.is_file():
        raise SystemExit(f"[ERROR] resolve_water: {path} missing — this FF defines no water models")
    names = [ln.split()[0] for ln in path.read_text().splitlines() if ln.strip() and not ln.lstrip().startswith(";")]
    if not names:
        raise SystemExit(f"[ERROR] resolve_water: {path} lists no water models")
    return names


def count_sites(itp: Path) -> int:
    """Atoms in the first [ atoms ] block of the water .itp (virtual sites count: they are atoms)."""
    if not itp.is_file():
        raise SystemExit(f"[ERROR] resolve_water: {itp} missing")
    in_atoms = False
    n = 0
    for ln in itp.read_text().splitlines():
        text = ln.split(";", 1)[0].strip()
        if not text or text.startswith("#"):
            continue
        m = re.fullmatch(r"\[\s*(\w+)\s*\]", text)
        if m:
            if in_atoms:
                break
            in_atoms = m.group(1).lower() == "atoms"
            continue
        if in_atoms:
            n += 1
    if n == 0:
        raise SystemExit(f"[ERROR] resolve_water: no [ atoms ] in {itp}")
    return n


def main(argv: list[str]) -> int:
    if len(argv) < 4:
        raise SystemExit(__doc__)
    ff, water, dirs = argv[1], argv[2], argv[3:]
    ffdir = find_ffdir(ff, dirs)
    menu = water_menu(ffdir)
    if water not in menu:
        raise SystemExit(f"[ERROR] resolve_water: WATER='{water}' is not a water model of {ffdir.name} "
                         f"(watermodels.dat lists: {' '.join(menu)})")
    if water in PDB2GMX_BUILTIN:
        arg, answer = water, "-"
    else:
        arg, answer = "select", str(menu.index(water) + 1)
    sites = count_sites(ffdir / f"{water}.itp")
    own_box = ffdir / f"{water}.gro"
    if own_box.is_file():
        box = str(own_box)
    elif sites in BOX_BY_SITES:
        box = BOX_BY_SITES[sites]
    else:
        raise SystemExit(f"[ERROR] resolve_water: no solvent box for a {sites}-site water ({water})")
    print(arg, answer, sites, box)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
