#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# plot_dssp.py — secondary-structure map (residue × frame) from gmx dssp output
# ─────────────────────────────────────────────────────────────────────────────
# Trajectory-agnostic: renders the .dat written by `gmx dssp -o` for either the
# plain-MD or T-REMD pipeline. The .dat has one row per frame, each row a string
# of single-letter DSSP codes (one per residue).
#
# Produces a categorical heatmap: residue index (y) vs frame (x), coloured by
# secondary-structure type, with a legend.
#
# Usage:
#   python plot_dssp.py DSSP_DAT OUT_PNG
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import sys
from pathlib import Path

# Human-readable names for the standard DSSP one-letter codes. Codes not in this
# map (e.g. version-specific additions) are still plotted, labelled by their
# raw letter — so this never crashes on an unexpected code.
DSSP_NAMES = {
    "H": "α-helix",
    "G": "3₁₀-helix",
    "I": "π-helix",
    "P": "PP-II helix",
    "E": "β-strand",
    "B": "β-bridge",
    "T": "turn",
    "S": "bend",
    "~": "coil",
    " ": "coil",
    "C": "coil",
    "-": "coil",
}

# Fixed, semantic colours grouped by structure family so the map reads at a
# glance and stays consistent across runs (unlike a discovery-order colormap):
#   helices  → blues/purple, strands → warm reds, turns/bends → greens,
#   coil     → light grey, so structured regions pop against the loopy default.
# Codes absent from this map fall back to a spare-colour cycle (never crashes).
DSSP_COLORS = {
    "H": "#2c5fa8",   # α-helix — deep blue
    "G": "#5b8dd6",   # 3₁₀-helix — mid blue
    "I": "#8fb3e8",   # π-helix — light blue
    "P": "#7b5ea7",   # PP-II — purple
    "E": "#c0392b",   # β-strand — red
    "B": "#e8954e",   # β-bridge — orange
    "T": "#3aa66f",   # turn — green
    "S": "#7fc99a",   # bend — light green
    "~": "#dcdcdc",   # coil
    " ": "#dcdcdc",   # coil
    "C": "#dcdcdc",   # coil
    "-": "#dcdcdc",   # coil
}
# Spare colours for any DSSP code not in DSSP_COLORS (kept distinct from above).
DSSP_FALLBACK_COLORS = ["#d4ac0d", "#16a085", "#cb4335", "#717d7e", "#884ea0"]


def read_dssp(path: Path) -> list[str]:
    """Return one SS-code string per frame, dropping blank/comment lines."""
    frames = [
        line.rstrip("\n")
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    ]
    if not frames:
        sys.exit(f"[ERROR] No frames found in {path}")
    nres = len(frames[0])
    # ASSUMES: every frame reports the same number of residues.
    assert all(len(f) == nres for f in frames), (
        f"inconsistent residue counts across frames in {path}"
    )
    return frames


def main():
    ap = argparse.ArgumentParser(
        description="Plot a gmx dssp .dat as a residue×frame secondary-structure map."
    )
    ap.add_argument("dat", type=Path, help="input gmx dssp .dat file")
    ap.add_argument("png", type=Path, help="output .png file")
    args = ap.parse_args()

    if not args.dat.is_file():
        sys.exit(f"[ERROR] dssp .dat not found: {args.dat}")

    try:
        import matplotlib
        matplotlib.use("Agg")          # headless: no display on compute nodes
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
        from matplotlib.patches import Patch
        import numpy as np
    except ImportError:
        sys.exit("[ERROR] matplotlib/numpy not available; cannot plot.")

    frames = read_dssp(args.dat)
    n_frames = len(frames)
    n_res = len(frames[0])

    # Discover the SS codes actually present and assign each an integer category.
    codes = sorted({c for frame in frames for c in frame})
    code_to_int = {c: i for i, c in enumerate(codes)}

    # grid[residue, frame] = category integer
    grid = np.array(
        [[code_to_int[frame[r]] for frame in frames] for r in range(n_res)]
    )

    # Discrete colormap: one fixed semantic colour per present code, falling back
    # to a spare cycle for any code outside DSSP_COLORS.
    fallback = iter(DSSP_FALLBACK_COLORS)
    color_list = [
        DSSP_COLORS.get(c) or next(fallback, "#000000") for c in codes
    ]
    cmap = ListedColormap(color_list)
    norm = BoundaryNorm(range(len(codes) + 1), cmap.N)

    # Bound the figure size: imshow(aspect="auto") stretches the data to fill the
    # axes, so the figure need not grow with frame count. Without the caps the width
    # scaled unbounded (e.g. 25k frames → a 500" sliver).
    fig_w = min(max(8.0, n_frames * 0.02), 16.0)
    fig_h = min(max(4.0, n_res * 0.05), 16.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(grid, aspect="auto", origin="lower", cmap=cmap, norm=norm,
              interpolation="nearest",
              extent=[0, n_frames, 1, n_res + 1])
    ax.set_xlabel("Frame")
    ax.set_ylabel("Residue")
    ax.set_title(f"Secondary structure ({n_res} residues, {n_frames} frames)")

    legend_handles = [
        Patch(facecolor=cmap(code_to_int[c]),
              label=f"{c if c.strip() else '·'}  {DSSP_NAMES.get(c, 'other')}")
        for c in codes
    ]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.01, 1.0),
              loc="upper left", fontsize=8, title="DSSP")

    fig.tight_layout()
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, dpi=150)
    plt.close(fig)
    print(f"[OK] Plot written to: {args.png}")


if __name__ == "__main__":
    main()
