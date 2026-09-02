"""Tests for the iterative NPT density-equilibration plateau test.

density_converged.py lives in scripts/simulation/ rather than in gromd_analysis
because it must be stdlib-only: the engines call it under the sbatch's system
python3, before the analysis conda env is activated. So it is loaded here by
path instead of imported as a package module.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "simulation" / "density_converged.py"


def _load():
    spec = importlib.util.spec_from_file_location("density_converged", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


density_converged = _load()


def _run(tol, min_seg, volumes):
    """Invoke the CLI the way the engines do; return (verdict:int, reason:str).

    One stdout line, "<verdict> <reason>" — the engines split it the same way, so
    the reason lands in the job log beside the per-segment volumes instead of in
    the separate .err stream.
    """
    proc = subprocess.run(
        [sys.executable, str(_SCRIPT), "--tol", str(tol), "--min-seg", str(min_seg),
         *(str(v) for v in volumes)],
        capture_output=True, text=True, check=True,
    )
    assert proc.stderr == "", "the reason belongs on stdout, not stderr"
    verdict, _, reason = proc.stdout.strip().partition(" ")
    return int(verdict), reason


def test_flat_series_converges():
    """Pure scatter about a constant mean has ~zero slope."""
    volumes = [500.0, 499.8, 500.2, 500.1, 499.9, 500.0, 500.3, 499.7]
    assert _run(0.005, 8, volumes)[0] == 1


def test_steady_drift_is_rejected_though_each_step_is_under_tol():
    """The regression that motivated replacing the consecutive-segment test.

    Each step here is 2/518 = 0.39 % — under the 0.5 % tolerance, so the old
    |V_n - V_{n-1}| / V_{n-1} test passed on every single comparison — while the
    box contracts 2.7 % across the window.
    """
    volumes = [520.0, 518.0, 516.0, 514.0, 512.0, 510.0, 508.0, 506.0]
    for prev, curr in zip(volumes, volumes[1:], strict=False):
        assert abs(curr - prev) / prev < 0.005, "premise: the old test would have passed"
    assert _run(0.005, 8, volumes)[0] == 0


def test_too_few_segments_is_not_converged():
    verdict, reason = _run(0.005, 8, [500.0, 499.0, 501.0])
    assert verdict == 0
    assert "too few" in reason


def test_reason_accompanies_the_verdict_on_one_stdout_line():
    verdict, reason = _run(0.005, 8, [520.0, 518.0, 516.0, 514.0, 512.0, 510.0, 508.0, 506.0])
    assert verdict == 0
    assert "still drifting" in reason and "2.7" in reason


def test_only_the_trailing_window_counts():
    """A large early transient must not keep a settled box from converging."""
    volumes = [600.0, 570.0, 540.0, 515.0] + [500.0, 500.1, 499.9, 500.0, 500.2, 499.8, 500.1, 499.9]
    assert _run(0.005, 8, volumes)[0] == 1


def test_direction_does_not_matter():
    """Expansion is rejected on the same footing as contraction."""
    volumes = [506.0, 508.0, 510.0, 512.0, 514.0, 516.0, 518.0, 520.0]
    assert _run(0.005, 8, volumes)[0] == 0


def test_drift_is_measured_across_the_whole_window():
    mod = density_converged
    # Exactly 1 nm^3 per segment over an 8-segment window about a mean of 500:
    # 7 nm^3 of drift / 500 = 1.4 %.
    volumes = [496.5, 497.5, 498.5, 499.5, 500.5, 501.5, 502.5, 503.5]
    test = mod.drift_over_window(volumes, 8)
    assert test.window == 8
    assert test.mean_volume == pytest.approx(500.0)
    assert test.drift_rel == pytest.approx(7.0 / 500.0)


def test_bad_arguments_exit_nonzero():
    for bad in (["--tol", "0", "--min-seg", "8"], ["--tol", "0.005", "--min-seg", "1"]):
        proc = subprocess.run([sys.executable, str(_SCRIPT), *bad, "500", "500"],
                              capture_output=True, text=True)
        assert proc.returncode != 0
