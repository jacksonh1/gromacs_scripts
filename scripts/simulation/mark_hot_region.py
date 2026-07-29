#!/usr/bin/env python3
"""Mark the REST2 solute (hot) region in a grompp -pp processed topology.

For every moleculetype whose name matches the solute pattern (default: names
starting with 'Protein'), append '_' to the atom-type column (column 2) of each
line in that moleculetype's [ atoms ] section. `plumed partial_tempering` reads
those '_' markers and scales the corresponding force-field terms.

Whole-protein, single-chain REST2 only: the solute is the protein, nothing else.
(Water and ions are left unmarked → unscaled, which is correct REST2.)

Reads the processed topology from stdin (or a file arg), writes the marked
topology to stdout. Fails loudly: exits 1 if zero atoms were marked.

Pure stdlib — runs under the sbatch's system python3.
"""

import re
import sys

SOLUTE_MOLTYPE = re.compile(r"^Protein")  # pdb2gmx names: Protein, Protein_chain_A, ...

_SECTION = re.compile(r"^\s*\[\s*(\S+)\s*\]")


def strip_comment(line):
    """Split a topology line into (data, comment) where comment includes the ';'."""
    i = line.find(";")
    if i == -1:
        return line, ""
    return line[:i], line[i:]


def mark(lines):
    out = []
    section = None          # current [ section ] name
    moltype = None          # current moleculetype name (None until named)
    awaiting_name = False    # just saw [ moleculetype ], next data line is "Name nrexcl"
    n_marked = 0

    for line in lines:
        raw = line.rstrip("\n")
        data, comment = strip_comment(raw)
        stripped = data.strip()

        m = _SECTION.match(data)
        if m:
            section = m.group(1)
            if section == "moleculetype":
                awaiting_name = True
                moltype = None
            out.append(raw)
            continue

        # Blank or comment-only line — passthrough, no state change.
        if stripped == "":
            out.append(raw)
            continue

        if awaiting_name:
            # First data line of a [ moleculetype ]: "<name> <nrexcl>"
            moltype = stripped.split()[0]
            awaiting_name = False
            out.append(raw)
            continue

        if section == "atoms" and moltype is not None and SOLUTE_MOLTYPE.match(moltype):
            fields = stripped.split()
            # [ atoms ] format: nr type resnr residue atom cgnr charge [mass]
            assert len(fields) >= 6, f"malformed [atoms] line in {moltype}: {raw!r}"
            fields[1] = fields[1] + "_"
            out.append("  ".join(fields) + ("  " + comment if comment else ""))
            n_marked += 1
            continue

        out.append(raw)

    if n_marked == 0:
        sys.stderr.write(
            "[ERROR] mark_hot_region: marked 0 atoms — no moleculetype matched "
            f"/{SOLUTE_MOLTYPE.pattern}/. Is this a grompp -pp processed topology?\n"
        )
        sys.exit(1)

    sys.stderr.write(f"[OK] mark_hot_region: marked {n_marked} solute atoms\n")
    return out


def main():
    src = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    lines = src.readlines()
    for line in mark(lines):
        sys.stdout.write(line + "\n")


if __name__ == "__main__":
    main()
