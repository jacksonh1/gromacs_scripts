"""Tests for topology parsing: [ molecules ], per-moltype atom counts, chain ranges."""

import pytest

from gromd_analysis.chains import (
    moltype_atom_counts,
    parse_molecules,
    protein_chains,
    read_gro_atom_names,
)

TOP = """\
; a system topology
#include "amber99sb-ildn.ff/forcefield.itp"

[ moleculetype ]
; name  nrexcl
Protein_chain_A   3

[ atoms ]
     1  N   1  PHE  N   1  -0.5   14.01
     2  CA  1  PHE  CA  2   0.0   12.01
     3  C   1  PHE  C   3   0.6   12.01

[ moleculetype ]
Protein_chain_B   3

[ atoms ]
     1  N   1  ALA  N   1  -0.5   14.01
     2  CA  1  ALA  CA  2   0.0   12.01

[ system ]
Two chains in water

[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     2      ; two copies of B
SOL              5000
NA                  7
"""


@pytest.fixture
def build_dir(tmp_path):
    (tmp_path / "system.top").write_text(TOP)
    return tmp_path


def test_parse_molecules_reads_only_the_molecules_section(build_dir):
    mols = parse_molecules(build_dir / "system.top")
    assert mols == [
        ("Protein_chain_A", 1),
        ("Protein_chain_B", 2),
        ("SOL", 5000),
        ("NA", 7),
    ]


def test_parse_molecules_strips_trailing_comments(build_dir):
    # 'Protein_chain_B 2 ; two copies of B' must parse as a count of 2, not fail.
    assert dict(parse_molecules(build_dir / "system.top"))["Protein_chain_B"] == 2


def test_moltype_atom_counts_counts_per_moleculetype(build_dir):
    counts = moltype_atom_counts(build_dir)
    assert counts["Protein_chain_A"] == 3
    assert counts["Protein_chain_B"] == 2


def test_protein_chains_expands_copies_and_assigns_contiguous_ranges(build_dir):
    chains = protein_chains(build_dir / "system.top", build_dir)
    # One copy of A (3 atoms) then two copies of B (2 atoms each) = 3 chains.
    assert [c.letter for c in chains] == ["A", "B", "C"]
    assert [(c.start, c.end) for c in chains] == [(1, 3), (4, 5), (6, 7)]


def test_protein_chains_ignores_solvent_and_ions(build_dir):
    assert len(protein_chains(build_dir / "system.top", build_dir)) == 3


def test_protein_chains_fails_loudly_on_unknown_moltype(tmp_path):
    # A [ molecules ] entry with no matching [ moleculetype ] must not be silently skipped.
    (tmp_path / "system.top").write_text(
        "[ molecules ]\nProtein_chain_Z  1\n"
    )
    with pytest.raises(SystemExit):
        protein_chains(tmp_path / "system.top", tmp_path)


def test_read_gro_atom_names(tmp_path):
    gro = tmp_path / "x.gro"
    gro.write_text(
        "a title\n"
        "    3\n"
        "    1PHE      N    1   0.000   0.000   0.000\n"
        "    1PHE     CA    2   0.100   0.000   0.000\n"
        "    1PHE      C    3   0.200   0.000   0.000\n"
        "   1.0   1.0   1.0\n"
    )
    assert read_gro_atom_names(gro) == ["N", "CA", "C"]
