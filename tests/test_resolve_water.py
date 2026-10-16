"""Tests for scripts/simulation/resolve_water.py (STEP 1 water-model resolution).

Invoked as a CLI by path, like the engines do. Force-field directories are synthesized in
tmp_path with just the files the resolver reads (watermodels.dat, <water>.itp, <water>.gro).
"""

import subprocess
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "simulation" / "resolve_water.py"

ITP_3SITE = "[ moleculetype ]\nSOL 2\n[ atoms ]\n1 OW 1 SOL OW 1 -0.834 16.0\n2 HW 1 SOL HW1 1 0.417 1.008\n3 HW 1 SOL HW2 1 0.417 1.008\n[ settles ]\n1 1 0.09572 0.15139\n"
ITP_4SITE = ITP_3SITE.replace("[ settles ]", "4 MW 1 SOL MW 1 0 0\n[ settles ]")


def make_ff(root, name, models, own_box=None):
    ff = root / f"{name}.ff"
    ff.mkdir(parents=True)
    (ff / "watermodels.dat").write_text("".join(f"{m}  {m.upper()}  some water\n" for m, _ in models))
    for m, itp in models:
        (ff / f"{m}.itp").write_text(itp)
    if own_box:
        (ff / f"{own_box}.gro").write_text("box\n0\n1 1 1\n")
    return ff


def run(ff, water, *dirs):
    return subprocess.run([sys.executable, str(_SCRIPT), ff, water, *map(str, dirs)],
                          capture_output=True, text=True)


def test_builtin_3site_is_unchanged(tmp_path):
    """tip3p must resolve to exactly the old engine behaviour: -water tip3p, spc216.gro."""
    make_ff(tmp_path, "amber", [("tip3p", ITP_3SITE), ("tip4pew", ITP_4SITE)])
    proc = run("amber", "tip3p", tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["tip3p", "-", "3", "spc216.gro"]


def test_builtin_4site_gets_4site_box(tmp_path):
    make_ff(tmp_path, "amber", [("tip3p", ITP_3SITE), ("tip4pew", ITP_4SITE)])
    assert run("amber", "tip4pew", tmp_path).stdout.split() == ["tip4pew", "-", "4", "tip4p.gro"]


def test_ff_specific_water_uses_select_and_own_box(tmp_path):
    """a99SB-disp style: not a pdb2gmx built-in -> `-water select` + menu index, FF's own box."""
    ff = make_ff(tmp_path, "disp", [("spc", ITP_3SITE), ("a99SBdisp_water", ITP_4SITE)],
                 own_box="a99SBdisp_water")
    out = run("disp", "a99SBdisp_water", tmp_path).stdout.split()
    assert out == ["select", "2", "4", str(ff / "a99SBdisp_water.gro")]


def test_search_order_first_dir_wins(tmp_path):
    make_ff(tmp_path / "a", "amber", [("tip3p", ITP_3SITE)])
    make_ff(tmp_path / "b", "amber", [("tip4p", ITP_4SITE)])
    assert run("amber", "tip3p", tmp_path / "a", tmp_path / "b").returncode == 0
    assert run("amber", "tip4p", tmp_path / "a", tmp_path / "b").returncode != 0


def test_water_not_in_ff_fails(tmp_path):
    make_ff(tmp_path, "disp", [("a99SBdisp_water", ITP_4SITE)])
    proc = run("disp", "tip3p", tmp_path)
    assert proc.returncode != 0
    assert "not a water model" in proc.stderr


def test_missing_ff_fails(tmp_path):
    proc = run("nosuch", "tip3p", tmp_path)
    assert proc.returncode != 0
    assert "not found" in proc.stderr


def test_missing_itp_fails(tmp_path):
    ff = make_ff(tmp_path, "amber", [("tip3p", ITP_3SITE)])
    (ff / "tip3p.itp").unlink()
    assert run("amber", "tip3p", tmp_path).returncode != 0
