"""Tests for the gmx dssp .dat reader."""

import pytest

from gromd_analysis.dssp import read_dssp


def test_reads_one_string_per_frame(tmp_path):
    p = tmp_path / "dssp.dat"
    p.write_text("# a comment\n@ directive\nHHHH-\nHHH--\n\nEEEE-\n")
    assert read_dssp(p) == ["HHHH-", "HHH--", "EEEE-"]


def test_inconsistent_residue_counts_are_rejected(tmp_path):
    p = tmp_path / "dssp.dat"
    p.write_text("HHHH-\nHHH\n")
    with pytest.raises(AssertionError):
        read_dssp(p)


def test_empty_file_exits(tmp_path):
    p = tmp_path / "dssp.dat"
    p.write_text("# only comments\n")
    with pytest.raises(SystemExit):
        read_dssp(p)
