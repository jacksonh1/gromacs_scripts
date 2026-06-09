#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# multichain_chain_index.py — build a protein-only index with per-chain groups
# ─────────────────────────────────────────────────────────────────────────────
# Owns all chain logic for the multi-chain analysis pipeline. From the system
# topology it determines the ordered protein chains and their atom counts, then —
# using a protein-only structure for atom names — writes a GROMACS index file with,
# for each chain:
#     [ ChainA ]            all atoms of the chain
#     [ ChainA_Backbone ]   its backbone atoms (N, CA, C, O)
#
# Chains are the protein MOLECULE TYPES from the topology (the physically-correct
# split), NOT `gmx splitch` — splitch splits on chain-ID/numbering jumps and can
# over-split a single moltype (observed 2 moltypes -> 3 "chains" on an internal gap).
#
# The atom indices are 1-based and match the protein-only trajectory's atom order
# (chains laid out in [ molecules ] order), so the same index serves the init
# reference and the stripped/aligned trajectory.
#
# Usage:
#   python3 multichain_chain_index.py BUILD_DIR PROTEIN_GRO OUT_NDX
#
#   BUILD_DIR    the job's build/ directory (holds <base>.top + chain *.itp)
#   PROTEIN_GRO  protein-only structure (e.g. <prefix>_stripped_aligned.gro)
#   OUT_NDX      index file to write
#
# Prints a final line "CHAINS: A B ..." (the chain letters) for the caller to parse.
# Exits non-zero (and prints nothing on that line) if <2 protein chains are found.
# ─────────────────────────────────────────────────────────────────────────────

import string
import sys
from dataclasses import dataclass
from pathlib import Path

BACKBONE_NAMES = {"N", "CA", "C", "O"}


@dataclass
class Chain:
    letter: str
    start: int          # 1-based, inclusive
    end: int            # 1-based, inclusive


def parse_molecules(top: Path) -> list[tuple[str, int]]:
    """Ordered (moltype, count) from the topology's [ molecules ] section."""
    mols: list[tuple[str, int]] = []
    in_section = False
    for raw in top.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            in_section = line.strip("[] ").lower() == "molecules"
            continue
        if in_section:
            parts = line.split()
            # ASSUMES: "<moltype> <count>" rows (gmx topology format)
            mols.append((parts[0], int(parts[1])))
    return mols


def moltype_atom_counts(build_dir: Path) -> dict[str, int]:
    """moltype name -> atom count, scanning every .top/.itp in build_dir.

    A file may hold several [ moleculetype ] blocks; we count [ atoms ] rows for
    whichever moleculetype is currently in scope. Files without a moleculetype
    (e.g. position-restraint itps) contribute nothing.
    """
    counts: dict[str, int] = {}
    for path in sorted(list(build_dir.glob("*.itp")) + list(build_dir.glob("*.top"))):
        section = None
        current: str | None = None
        for raw in path.read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                if section == "moleculetype":
                    current = None
                continue
            if section == "moleculetype" and current is None:
                current = line.split()[0]
            elif section == "atoms" and current is not None:
                counts[current] = counts.get(current, 0) + 1
    return counts


def protein_chains(top: Path, build_dir: Path) -> list[Chain]:
    """Ordered protein chains with 1-based atom ranges in the protein-only system."""
    mols = parse_molecules(top)
    atom_counts = moltype_atom_counts(build_dir)
    chains: list[Chain] = []
    cursor = 1
    for moltype, count in mols:
        if not moltype.startswith("Protein"):
            continue
        if moltype not in atom_counts:
            sys.exit(f"[ERROR] No atom count found for protein moltype '{moltype}' "
                     f"in {build_dir} (.itp/.top). Cannot build chain index.")
        n = atom_counts[moltype]
        for _ in range(count):
            letter = string.ascii_uppercase[len(chains)]
            chains.append(Chain(letter, cursor, cursor + n - 1))
            cursor += n
    return chains


def read_gro_atom_names(gro: Path) -> list[str]:
    """Atom names (stripped) for each atom in a .gro, in file order."""
    lines = gro.read_text().splitlines()
    n_atoms = int(lines[1].split()[0])
    names = []
    for line in lines[2:2 + n_atoms]:
        # gro fixed columns: resnr(5) resname(5) atomname(5) atomnr(5) x y z
        names.append(line[10:15].strip())
    return names


def write_ndx(path: Path, groups: dict[str, list[int]]) -> None:
    with path.open("w") as fh:
        for name, idx in groups.items():
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(f"{a:5d}" for a in idx[i:i + 15]) + "\n")
            fh.write("\n")


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit("Usage: python3 multichain_chain_index.py BUILD_DIR PROTEIN_GRO OUT_NDX")
    build_dir, protein_gro, out_ndx = (Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))

    tops = list(build_dir.glob("*.top"))
    if not tops:
        sys.exit(f"[ERROR] No .top found in {build_dir}")
    chains = protein_chains(tops[0], build_dir)

    if len(chains) < 2:
        sys.exit(f"[ERROR] Found {len(chains)} protein chain(s); the multi-chain "
                 f"index needs >= 2. (Single-chain jobs use the standard path.)")

    names = read_gro_atom_names(protein_gro)
    n_expected = chains[-1].end
    # Fail loudly if the topology chain ranges don't match the protein-only structure.
    assert len(names) == n_expected, (
        f"protein structure has {len(names)} atoms but topology chains sum to "
        f"{n_expected} — atom ordering mismatch, refusing to build a wrong index"
    )

    groups: dict[str, list[int]] = {}
    for ch in chains:
        all_idx = list(range(ch.start, ch.end + 1))
        bb_idx = [i for i in all_idx if names[i - 1] in BACKBONE_NAMES]
        groups[f"Chain{ch.letter}"] = all_idx
        groups[f"Chain{ch.letter}_Backbone"] = bb_idx

    write_ndx(out_ndx, groups)

    print(f"[OK] Wrote chain index: {out_ndx}")
    for ch in chains:
        n_bb = len(groups[f'Chain{ch.letter}_Backbone'])
        print(f"     Chain{ch.letter}: atoms {ch.start}-{ch.end} "
              f"({ch.end - ch.start + 1} atoms, {n_bb} backbone)")
    # Machine-readable summary line for the caller.
    print("CHAINS: " + " ".join(ch.letter for ch in chains))


if __name__ == "__main__":
    main()
