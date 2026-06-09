#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# plot_xvg.py — generic line plot for a GROMACS .xvg file
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: reads the title and axis labels embedded in the .xvg
# itself (the `@ title` / `@ xaxis label` / `@ yaxis label` grace directives),
# so the same script renders RMSD, Rg, RMSF, and any other xvg from either the
# plain-MD or T-REMD pipeline with no changes.
#
# Plots every data column (>= column 1) against column 0, using the per-series
# legends from the xvg (`@ sN legend ...`) when present.
#
# Usage:
#   python plot_xvg.py PLOT_XVG OUT_PNG
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import re
import sys
from pathlib import Path
from dataclasses import dataclass, field


@dataclass
class XvgData:
    """Parsed contents of a GROMACS .xvg file."""
    title: str
    xlabel: str
    ylabel: str
    legends: dict[int, str]          # series index (0-based) -> legend text
    columns: list[list[float]]       # columns[0] = x, columns[1:] = y series


def _nm_to_angstrom(label: str) -> tuple[str, bool]:
    """GROMACS reports distances in nm; structural biologists read ångström.

    If an axis label carries an `nm` unit, relabel it to Å and signal that the
    corresponding data column must be scaled ×10. Returns (new_label, scale?).
    """
    if re.search(r"\bnm\b", label):
        return re.sub(r"\bnm\b", "Å", label), True
    return label, False


def parse_xvg(path: Path) -> XvgData:
    title = ""
    xlabel = ""
    ylabel = ""
    legends: dict[int, str] = {}
    rows: list[list[float]] = []

    title_re = re.compile(r'@\s+title\s+"(.*)"')
    xaxis_re = re.compile(r'@\s+xaxis\s+label\s+"(.*)"')
    yaxis_re = re.compile(r'@\s+yaxis\s+label\s+"(.*)"')
    legend_re = re.compile(r'@\s+s(\d+)\s+legend\s+"(.*)"')

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            if (m := title_re.match(line)):
                title = m.group(1)
            elif (m := xaxis_re.match(line)):
                xlabel = m.group(1)
            elif (m := yaxis_re.match(line)):
                ylabel = m.group(1)
            elif (m := legend_re.match(line)):
                legends[int(m.group(1))] = m.group(2)
            continue
        # Data line.
        rows.append([float(tok) for tok in line.split()])

    if not rows:
        sys.exit(f"[ERROR] No data rows found in {path}")

    ncol = len(rows[0])
    # ASSUMES: every data row has the same column count (gmx xvg output does).
    assert all(len(r) == ncol for r in rows), f"ragged columns in {path}"
    columns = [[r[c] for r in rows] for c in range(ncol)]

    return XvgData(title, xlabel, ylabel, legends, columns)


def main():
    ap = argparse.ArgumentParser(description="Plot a GROMACS .xvg as a line chart.")
    ap.add_argument("xvg", type=Path, help="input .xvg file")
    ap.add_argument("png", type=Path, help="output .png file")
    args = ap.parse_args()

    if not args.xvg.is_file():
        sys.exit(f"[ERROR] xvg not found: {args.xvg}")

    try:
        import matplotlib
        matplotlib.use("Agg")          # headless: no display on compute nodes
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("[ERROR] matplotlib not available; cannot plot. "
                 "Install matplotlib or run without plotting.")

    data = parse_xvg(args.xvg)

    # Convert any nm axis to ångström so RMSD/Rg/RMSF plot in the units expected.
    data.xlabel, x_is_nm = _nm_to_angstrom(data.xlabel)
    data.ylabel, y_is_nm = _nm_to_angstrom(data.ylabel)
    if x_is_nm:
        data.columns[0] = [v * 10.0 for v in data.columns[0]]
    if y_is_nm:
        data.columns[1:] = [[v * 10.0 for v in col] for col in data.columns[1:]]

    x = data.columns[0]
    y_series = data.columns[1:]
    if not y_series:
        sys.exit(f"[ERROR] {args.xvg} has only one column; nothing to plot.")

    fig, ax = plt.subplots(figsize=(8, 4))
    for i, y in enumerate(y_series):
        # legends are keyed by series index (0-based), matching y_series order
        label = data.legends.get(i)
        ax.plot(x, y, linewidth=0.8, label=label)

    ax.set_xlabel(data.xlabel or "x")
    ax.set_ylabel(data.ylabel or "y")
    ax.set_title(data.title or args.xvg.stem)
    # Only show a legend if at least one series was labelled.
    if any(data.legends.get(i) for i in range(len(y_series))):
        ax.legend(fontsize=8)
    fig.tight_layout()
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, dpi=150)
    plt.close(fig)
    print(f"[OK] Plot written to: {args.png}")


if __name__ == "__main__":
    main()
