"""Tests for the STEP 7 REST2 topology check (scripts/simulation/check_rest2_topology.py).

The checker runs under the sbatch's system python3 like the other STEP helpers, so it is
invoked here as a CLI by path. Topologies are synthesized to look like
`plumed partial_tempering` output (hot types carry the `_` suffix, dihedrals inline);
each failure test perturbs one thing a broken scaling would get wrong.
"""

import math
import subprocess
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "simulation" / "check_rest2_topology.py"


def make_top(lam, *, charge_factor=None, dih_factor=None, nb_factor=None, extra="", dih_funct=9,
             bond_b0="0.1010"):
    """A two-atom-type solute (hot) + one water type (cold), scaled to `lam` like partial_tempering."""
    sq = math.sqrt(lam)
    cf = sq if charge_factor is None else charge_factor
    df = lam if dih_factor is None else dih_factor
    nf = sq if nb_factor is None else nb_factor
    dih_params = "0.0 {k:.6g} 3" if dih_funct == 9 else "{k:.6g}"
    k = 0.65084 * df
    return f"""\
[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
OB 8 16.00 0.0000 A 0.295992 0.87864
OB_ OB 16.00 0.0000 A 0.295992 {0.87864 * lam:.6g}  ; scaled
HB 1 1.008 0.0000 A 0.106908 0.0656888
HB_ HB 1.008 0.0000 A 0.106908 {0.0656888 * lam:.6g}  ; scaled
OW 8 16.00 0.0000 A 0.3165 0.998989
OW_ OW 16.00 0.0000 A 0.3165 {0.998989 * lam:.6g}  ; scaled
[ nonbond_params ]
OB HB 1 0.150 1.2552
OB_ HB 1 0.15 {1.2552 * nf:.6g}  ; scaled
OB HB_ 1 0.15 {1.2552 * nf:.6g}  ; scaled
OB_ HB_ 1 0.15 {1.2552 * lam:.6g}  ; scaled
[ moleculetype ]
Protein 3
[ atoms ]
1 HB_ 1 ALA H 1 {0.2719 * cf:.6g} 1.008
2 OB_ 1 ALA O 2 {-0.5679 * cf:.6g} 16.00
3 HB_ 1 ALA H2 3 {0.2719 * cf:.6g} 1.008
4 OB_ 1 ALA O2 4 {-0.5679 * cf:.6g} 16.00
[ bonds ]
1 2 1 {bond_b0} 363171.2
[ pairs ]
1 4 1
[ dihedrals ]
1 2 3 4 {dih_funct} {dih_params.format(k=k)}
{extra}[ moleculetype ]
SOL 1
[ atoms ]
1 OW 1 SOL OW 1 0 16.00
[ system ]
test
[ molecules ]
Protein 1
SOL 10
"""


def _run(tmp_path, cold, hot, lam_hot=0.5):
    (tmp_path / "cold.top").write_text(cold)
    (tmp_path / "hot.top").write_text(hot)
    return subprocess.run(
        [sys.executable, str(_SCRIPT), str(tmp_path / "cold.top"), str(tmp_path / "hot.top"),
         "--lambda-cold", "1.0", "--lambda-hot", str(lam_hot)],
        capture_output=True, text=True,
    )


def test_correct_scaling_passes(tmp_path):
    proc = _run(tmp_path, make_top(1.0), make_top(0.5))
    assert proc.returncode == 0, proc.stderr
    assert "[OK]" in proc.stdout


def test_unscaled_charges_fail(tmp_path):
    proc = _run(tmp_path, make_top(1.0), make_top(0.5, charge_factor=1.0))
    assert proc.returncode == 1
    assert "charge" in proc.stderr


def test_unscaled_dihedral_fails(tmp_path):
    proc = _run(tmp_path, make_top(1.0), make_top(0.5, dih_factor=1.0))
    assert proc.returncode == 1
    assert "dihedral parameter" in proc.stderr


def test_wrong_nonbond_params_fails(tmp_path):
    """The a99SB-disp backbone O-H override: a mixed hot/cold pair must get eps*sqrt(lambda)."""
    proc = _run(tmp_path, make_top(1.0), make_top(0.5, nb_factor=0.5))
    assert proc.returncode == 1
    assert "nonbond_params" in proc.stderr


def test_cmap_is_rejected(tmp_path):
    cmap = "[ cmap ]\n1 2 3 4 1 1\n"
    proc = _run(tmp_path, make_top(1.0, extra=cmap), make_top(0.5, extra=cmap))
    assert proc.returncode == 1
    assert "[ cmap ]" in proc.stderr


def test_unscalable_dihedral_function_is_rejected(tmp_path):
    """funct 8 (tabulated) passes through partial_tempering unscaled."""
    proc = _run(tmp_path, make_top(1.0, dih_funct=8), make_top(0.5, dih_funct=8))
    assert proc.returncode == 1
    assert "funct 8" in proc.stderr


def test_changed_bond_fails(tmp_path):
    """REST2 must leave bonds alone; any difference means the tops don't match."""
    proc = _run(tmp_path, make_top(1.0), make_top(0.5, bond_b0="0.1020"))
    assert proc.returncode == 1
    assert "bonds" in proc.stderr


def test_unmarked_topology_fails(tmp_path):
    """No hot atoms at all (mark_hot_region skipped) must not pass as 'nothing to check'."""
    cold = make_top(1.0).replace("HB_ 1 ALA", "HB 1 ALA").replace("OB_ 1 ALA", "OB 1 ALA")
    proc = _run(tmp_path, cold, cold)
    assert proc.returncode == 1


def test_identical_lambda_ladder_passes(tmp_path):
    """The identical-lambda exchange sanity run (all replicas at lambda=1) is legitimate."""
    proc = _run(tmp_path, make_top(1.0), make_top(1.0), lam_hot=1.0)
    assert proc.returncode == 0, proc.stderr


def test_hot_file_scaled_but_lambda_claimed_one_fails(tmp_path):
    proc = _run(tmp_path, make_top(1.0), make_top(0.5), lam_hot=1.0)
    assert proc.returncode == 1


def test_lambda_order_is_validated(tmp_path):
    (tmp_path / "cold.top").write_text(make_top(1.0))
    (tmp_path / "hot.top").write_text(make_top(0.5))
    proc = subprocess.run([sys.executable, str(_SCRIPT), str(tmp_path / "cold.top"), str(tmp_path / "hot.top"),
                           "--lambda-cold", "0.5", "--lambda-hot", "1.0"], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "lambda-hot <= lambda-cold" in proc.stderr
