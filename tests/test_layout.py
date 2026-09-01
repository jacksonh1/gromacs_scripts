"""Tests for JobDir.detect — the parse-don't-validate boundary at the job directory."""

import pytest

from gromd_analysis.layout import JobDir, JobLayoutError

TOP_TWO_CHAINS = """\
[ moleculetype ]
Protein_chain_A  3
[ atoms ]
     1  N   1  PHE  N   1  -0.5   14.01
[ moleculetype ]
Protein_chain_B  3
[ atoms ]
     1  N   1  ALA  N   1  -0.5   14.01
[ molecules ]
Protein_chain_A  1
Protein_chain_B  1
SOL           3000
"""

TOP_ONE_CHAIN = """\
[ molecules ]
Protein_chain_A  1
SOL           3000
"""


def make_job(root, mode, rep="000", top=TOP_ONE_CHAIN, legacy_traj=False):
    """Build a minimal on-disk job tree. Contents don't matter — detect() checks existence."""
    stage = root / "prod" if mode == "MD" else root / "prod" / f"rep{rep}"
    stage.mkdir(parents=True)
    basename = {"MD": "md", "REMD": "remd", "REST2": "rest2"}[mode]
    (stage / f"{basename}.tpr").write_text("")
    if legacy_traj:
        (root / "trajectories").mkdir()
        (root / "trajectories" / f"{basename}.xtc").write_text("")
    else:
        (stage / f"{basename}.xtc").write_text("")
    (root / "em").mkdir()
    if top is not None:
        (root / "build").mkdir()
        (root / "build" / "system.top").write_text(top)
    return root


@pytest.mark.parametrize(
    "mode,basename,prefix_stem",
    [("MD", "md", "md"), ("REMD", "remd", "remd_rep000"), ("REST2", "rest2", "rest2_rep000")],
)
def test_detects_each_engine(tmp_path, mode, basename, prefix_stem):
    job = JobDir.detect(make_job(tmp_path, mode))
    assert job.mode == mode
    assert job.tpr.name == f"{basename}.tpr"
    assert job.xtc.name == f"{basename}.xtc"
    assert job.prefix.name == prefix_stem
    assert job.prefix.parent == job.analysis_dir


def test_replica_slot_is_honoured(tmp_path):
    job = JobDir.detect(make_job(tmp_path, "REMD", rep="007"), rep="007")
    assert job.tpr.parent.name == "rep007"
    assert job.prefix.name == "remd_rep007"


def test_detect_does_not_create_the_analysis_dir(tmp_path):
    # detect() is a parser: it resolves paths, it does not have side effects.
    job = JobDir.detect(make_job(tmp_path, "MD"))
    assert not job.analysis_dir.exists()


# ── failure modes: these must raise, never return a half-resolved JobDir ─────

def test_missing_directory_raises(tmp_path):
    with pytest.raises(JobLayoutError, match="Not a directory"):
        JobDir.detect(tmp_path / "nope")


def test_unfinished_job_raises(tmp_path):
    (tmp_path / "prod").mkdir()
    with pytest.raises(JobLayoutError, match="Cannot find"):
        JobDir.detect(tmp_path)


def test_purged_trajectory_raises(tmp_path):
    root = make_job(tmp_path, "REMD")
    (root / "prod" / "rep000" / "remd.xtc").unlink()
    with pytest.raises(JobLayoutError, match="Trajectory not found"):
        JobDir.detect(root)


def test_wrong_replica_slot_raises(tmp_path):
    with pytest.raises(JobLayoutError):
        JobDir.detect(make_job(tmp_path, "REMD"), rep="042")


# ── legacy layout ────────────────────────────────────────────────────────────

def test_legacy_trajectories_dir_is_accepted(tmp_path):
    job = JobDir.detect(make_job(tmp_path, "MD", legacy_traj=True))
    assert job.xtc.parent.name == "trajectories"


def test_legacy_fallback_is_not_used_when_the_stage_file_exists(tmp_path):
    root = make_job(tmp_path, "MD")
    (root / "trajectories").mkdir()
    (root / "trajectories" / "md.xtc").write_text("")
    assert JobDir.detect(root).xtc.parent.name == "prod"


# ── chain count drives the multichain dispatch ───────────────────────────────

def test_single_chain_topology(tmp_path):
    job = JobDir.detect(make_job(tmp_path, "MD"))
    assert job.n_chains == 1
    assert job.multichain is False


def test_multi_chain_topology(tmp_path):
    job = JobDir.detect(make_job(tmp_path, "MD", top=TOP_TWO_CHAINS))
    assert job.n_chains == 2
    assert job.multichain is True


def test_repeated_moltype_counts_every_copy(tmp_path):
    top = "[ molecules ]\nProtein_chain_A  3\nSOL  100\n"
    assert JobDir.detect(make_job(tmp_path, "MD", top=top)).n_chains == 3


def test_missing_topology_falls_back_to_one_chain_with_a_warning(tmp_path, capsys):
    job = JobDir.detect(make_job(tmp_path, "MD", top=None))
    assert job.n_chains == 1
    assert "[WARN]" in capsys.readouterr().err


def test_topology_without_protein_rows_warns(tmp_path, capsys):
    job = JobDir.detect(make_job(tmp_path, "MD", top="[ molecules ]\nSOL  100\n"))
    assert job.n_chains == 1
    assert "[WARN]" in capsys.readouterr().err


def test_chatter_never_goes_to_stdout(tmp_path, capsys):
    # gromd-layout's stdout is eval'd by run_analysis.sh — anything on it would be
    # executed as shell. Informational output must go to stderr.
    JobDir.detect(make_job(tmp_path, "MD", top=None, legacy_traj=True))
    assert capsys.readouterr().out == ""
