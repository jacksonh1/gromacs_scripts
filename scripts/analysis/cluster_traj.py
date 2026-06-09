#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
# cluster_traj.py — conformational clustering of a stripped/aligned trajectory
# ─────────────────────────────────────────────────────────────────────────────
# Groups the sampled frames into discrete conformational states (Cα-RMSD
# clustering) and reports populations + a representative structure per state.
# Trajectory-agnostic: serves both the plain-MD and T-REMD pipelines (it consumes
# the protein-only, PBC-fixed, backbone-aligned outputs of strip_and_align).
#
# Usage:
#   python cluster_traj.py STRUCT TRAJ OUT_PREFIX [options]
#
#   STRUCT  protein-only reference structure (.gro/.pdb); topology
#   TRAJ    protein-only trajectory (.xtc), already PBC-fixed and aligned
#   OUT_PREFIX  e.g. analysis/md or analysis/remd_rep000. Outputs are written to a
#               clustering/ subdir next to it, e.g.
#               analysis/clustering/remd_rep000_cluster_*
#
# Options:
#   --method {dbscan,kmeans}   default dbscan
#   --cutoff NM                RMSD cutoff in nm (default 0.20 = 2.0 Å); DBSCAN only
#   --min-samples N            DBSCAN min_samples; default auto = max(10, 1.5% of frames).
#                              The lever for cluster count vs noise — raise it for fewer,
#                              denser clusters (sparser regions fall into noise).
#   --n-clusters K             k-means cluster count (default 5)
#   --selection SEL            MDAnalysis selection to cluster on (default 'name CA')
#   --stride N                 use every Nth frame (default 1 = all frames)
#
# ── Why this and not `gmx cluster` ───────────────────────────────────────────
# `gmx cluster` (gromos) builds the full pairwise-RMSD matrix → O(N²) in time and
# memory, impractical for the long production runs (25k+ frames). This clusters on
# flattened Cα coordinates with scikit-learn, which scales to all frames: k-means
# is O(N); DBSCAN is heavier but still far cheaper than the O(N²) matrix.
#
# ── eps ↔ RMSD cutoff (DBSCAN) ───────────────────────────────────────────────
# The input frames are already aligned to one common backbone reference (by
# strip_and_align_trajectory.sh), so the Euclidean distance between two frames'
# flattened coordinate vectors equals √N_atoms × RMSD(relative to that common
# alignment). We therefore set the DBSCAN neighbourhood radius
#       eps = cutoff_Å × √N_selected_atoms
# so the user-facing knob is a genuine RMSD cutoff (in nm). This is the one fix
# carried over from the Amber cluster_MD.py, whose eps was a raw flattened-coord
# distance mislabelled as Å. ASSUMES: TRAJ is pre-aligned (it is, from the
# pipeline) — this script does NOT re-fit.
#
# ── Scaling caveat ───────────────────────────────────────────────────────────
# DBSCAN's neighbour graph grows large when nearly all frames fall within the
# cutoff (a very stable structure at a loose cutoff) — the trivial "one cluster"
# answer, but memory-heavy at tens of thousands of frames. Escape hatches:
# --stride, or --method kmeans (O(N), always scales).
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless: no display on compute nodes
import matplotlib.pyplot as plt

import MDAnalysis as mda
from sklearn.cluster import DBSCAN, KMeans

NOISE = -1   # DBSCAN label for frames in no dense region

# `min_samples` (how many neighbours within eps a frame needs to seed/join a cluster)
# is the DBSCAN density knob that decides what is a real state vs. noise. A fixed value
# does not scale — 10 frames is 5% of a 200-frame run but 0.1% of a 10k-frame run, which
# lets tiny knots become their own clusters (e.g. 80 clusters, most <1%). So the default
# scales with the trajectory length: a region must hold ~1.5% of the frames to be a
# cluster; anything sparser falls into noise. The 1.5% value was tuned on the WW-domain
# REMD slots (cutoff 0.20 nm) to keep even the most heterogeneous case to ≲10 clusters;
# raise --min-samples for even fewer.
MIN_SAMPLES_FRAC = 0.015   # default min_samples as a fraction of the clustered frames
MIN_SAMPLES_FLOOR = 10     # but never below this (keeps short runs sane)


def load_coords(struct: Path, traj: Path, selection: str, stride: int):
    """Return (coords, frame_indices, times_ps, universe).

    coords        : (n_used, n_sel*3) float64 — flattened selected-atom coords
    frame_indices : (n_used,) int — trajectory frame number of each row
    times_ps      : (n_used,) float — simulation time (ps) of each row
    """
    u = mda.Universe(str(struct), str(traj))
    atoms = u.select_atoms(selection)
    if atoms.n_atoms == 0:
        sys.exit(f"[ERROR] selection '{selection}' matched no atoms in {struct}")
    print(f"[INFO] clustering on {atoms.n_atoms} atoms (selection: '{selection}')")
    print(f"[INFO] {len(u.trajectory)} frames total; stride {stride}")

    coords, frame_indices, times_ps = [], [], []
    for ts in u.trajectory[::stride]:
        coords.append(atoms.positions.flatten())
        frame_indices.append(ts.frame)
        times_ps.append(float(ts.time))
    coords = np.asarray(coords, dtype=np.float64)
    if len(coords) < 2:
        sys.exit(f"[ERROR] need ≥2 frames to cluster, got {len(coords)}")
    return coords, np.asarray(frame_indices), np.asarray(times_ps), u


def cluster_dbscan(coords, cutoff_nm, min_samples):
    """DBSCAN with eps derived from the RMSD cutoff (see header). Returns labels."""
    n_atoms = coords.shape[1] // 3
    cutoff_ang = cutoff_nm * 10.0
    eps = cutoff_ang * np.sqrt(n_atoms)
    print(f"[INFO] DBSCAN: RMSD cutoff {cutoff_nm} nm ({cutoff_ang:.1f} Å) → "
          f"eps {eps:.2f} (= {cutoff_ang:.1f} Å × √{n_atoms}), min_samples {min_samples}")
    labels = DBSCAN(eps=eps, min_samples=min_samples,
                    metric="euclidean", n_jobs=-1).fit_predict(coords)
    n_clusters = len({l for l in labels if l != NOISE})
    if n_clusters == 0:
        sys.exit("[ERROR] DBSCAN found 0 clusters (all frames are noise). The trajectory "
                 "has no region dense enough to seed a cluster — lower --min-samples or "
                 "raise --cutoff.")
    return labels


def cluster_kmeans(coords, n_clusters):
    """k-means on flattened coords. Returns labels (no noise)."""
    print(f"[INFO] k-means: k={n_clusters}")
    return KMeans(n_clusters=n_clusters, random_state=42,
                  n_init=10).fit_predict(coords)


def rank_by_population(labels):
    """Relabel non-noise clusters 0,1,… by descending population (0 = largest).

    DBSCAN/k-means label numbers are arbitrary; ranking makes c00 the dominant
    state. Noise (-1) is preserved as -1.
    """
    real = labels[labels != NOISE]
    if real.size == 0:
        sys.exit("[ERROR] no non-noise frames to rank")
    uniq, counts = np.unique(real, return_counts=True)
    order = uniq[np.argsort(-counts)]                      # largest first
    remap = {old: new for new, old in enumerate(order)}
    remap[NOISE] = NOISE
    return np.array([remap[l] for l in labels])


def representative_frame(coords, labels, cluster_id):
    """Index (into the strided arrays) of the frame nearest the cluster centroid."""
    mask = labels == cluster_id
    members = coords[mask]
    centroid = members.mean(axis=0)
    local = int(np.argmin(np.linalg.norm(members - centroid, axis=1)))
    return int(np.where(mask)[0][local])


def write_assignments(out, frame_indices, times_ps, labels):
    np.savetxt(out, np.column_stack([frame_indices, times_ps, labels]),
               delimiter=",", header="frame,time_ps,cluster",
               fmt=["%d", "%.3f", "%d"], comments="")
    print(f"[OK] assignments → {out}")


def write_rep_pdbs(u, coords, labels, frame_indices, prefix):
    """Write the full-protein representative structure for each non-noise cluster."""
    reps = {}
    for cid in sorted(l for l in set(labels) if l != NOISE):
        idx = representative_frame(coords, labels, cid)
        u.trajectory[int(frame_indices[idx])]
        out = f"{prefix}_cluster_rep_c{cid:02d}.pdb"
        # A .gro topology carries no PDB metadata (chainIDs, occupancies, …); the
        # writer warns once per missing attr. Those are expected and cosmetic —
        # silence just this write so genuine warnings elsewhere still surface.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            u.atoms.write(out)
        reps[cid] = idx
        print(f"[OK] cluster c{cid:02d} rep → frame {frame_indices[idx]} → {out}")
    return reps


def write_summary(out, method, n_frames, labels, frame_indices, times_ps, reps,
                  cutoff_nm, min_samples, n_kmeans):
    lines = ["# Conformational clustering summary", f"method: {method}"]
    if method == "dbscan":
        lines += [f"rmsd_cutoff_nm: {cutoff_nm}", f"min_samples: {min_samples}"]
    else:
        lines += [f"n_clusters_requested: {n_kmeans}"]
    n_noise = int((labels == NOISE).sum())
    n_clusters = len({l for l in labels if l != NOISE})
    lines += [
        f"frames_clustered: {n_frames}",
        f"n_clusters: {n_clusters}",
        f"noise_frames: {n_noise} ({100*n_noise/n_frames:.1f}%)",
        "",
        f"{'cluster':>8} {'frames':>8} {'percent':>8} {'rep_frame':>10} {'rep_time_ps':>12}",
    ]
    for cid in sorted(l for l in set(labels) if l != NOISE):
        n = int((labels == cid).sum())
        idx = reps[cid]
        lines.append(f"{'c%02d' % cid:>8} {n:>8} {100*n/n_frames:>7.1f}% "
                     f"{frame_indices[idx]:>10} {times_ps[idx]:>12.1f}")
    if n_noise:
        lines.append(f"{'noise':>8} {n_noise:>8} {100*n_noise/n_frames:>7.1f}%")
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"[OK] summary → {out}")


def plot_populations(out, labels, n_frames):
    ids = sorted(set(labels), key=lambda c: (c == NOISE, c))   # clusters then noise
    counts = [int((labels == c).sum()) for c in ids]
    names = ["noise" if c == NOISE else f"c{c:02d}" for c in ids]
    colors = ["lightcoral" if c == NOISE else "steelblue" for c in ids]
    fig, ax = plt.subplots(figsize=(max(6, len(ids)), 5))
    ax.bar(names, counts, color=colors, edgecolor="black")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Frames")
    ax.set_title("Cluster populations")
    for i, c in enumerate(counts):
        ax.text(i, c, f"{100*c/n_frames:.1f}%", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[OK] populations plot → {out}")


def plot_timeseries(out, times_ps, labels):
    t_ns = times_ps / 1000.0
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(t_ns, labels, s=3, c=labels, cmap="tab10")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Cluster (−1 = noise)")
    ax.set_title("Cluster assignment over time")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[OK] timeseries plot → {out}")


def main():
    ap = argparse.ArgumentParser(
        description="Conformational clustering of a stripped/aligned trajectory.")
    ap.add_argument("struct", type=Path, help="protein-only topology (.gro/.pdb)")
    ap.add_argument("traj", type=Path, help="protein-only aligned trajectory (.xtc)")
    ap.add_argument("prefix", type=str,
                    help="output prefix; files go in <prefix-dir>/clustering/"
                         "<prefix-name>_cluster_*")
    ap.add_argument("--method", choices=["dbscan", "kmeans"], default="dbscan")
    ap.add_argument("--cutoff", type=float, default=0.20,
                    help="RMSD cutoff in nm (DBSCAN; default 0.20 = 2.0 Å)")
    ap.add_argument("--min-samples", type=int, default=None,
                    help="DBSCAN min_samples (frames a cluster's core point needs within "
                         "the cutoff). Default: auto = max(%d, %g%% of frames); raise it "
                         "for fewer clusters / more noise."
                         % (MIN_SAMPLES_FLOOR, MIN_SAMPLES_FRAC * 100))
    ap.add_argument("--n-clusters", type=int, default=5, help="k-means cluster count")
    ap.add_argument("--selection", type=str, default="name CA",
                    help="MDAnalysis selection to cluster on (default 'name CA')")
    ap.add_argument("--stride", type=int, default=1, help="use every Nth frame")
    args = ap.parse_args()

    if not args.struct.is_file():
        sys.exit(f"[ERROR] structure not found: {args.struct}")
    if not args.traj.is_file():
        sys.exit(f"[ERROR] trajectory not found: {args.traj}")

    coords, frame_indices, times_ps, u = load_coords(
        args.struct, args.traj, args.selection, args.stride)

    # Resolve min_samples: auto-scale with the frame count unless given explicitly.
    n_frames = len(coords)
    if args.min_samples is None:
        min_samples = max(MIN_SAMPLES_FLOOR, round(MIN_SAMPLES_FRAC * n_frames))
        print(f"[INFO] min_samples auto = {min_samples} "
              f"({MIN_SAMPLES_FRAC*100:g}% of {n_frames} frames, floor {MIN_SAMPLES_FLOOR})")
    else:
        min_samples = args.min_samples
        print(f"[INFO] min_samples = {min_samples} (user-set)")

    if args.method == "dbscan":
        labels = cluster_dbscan(coords, args.cutoff, min_samples)
    else:
        labels = cluster_kmeans(coords, args.n_clusters)
    labels = rank_by_population(labels)

    n_clusters = len({l for l in labels if l != NOISE})
    n_noise = int((labels == NOISE).sum())
    print(f"[INFO] {n_clusters} clusters, {n_noise} noise frames "
          f"({100*n_noise/n_frames:.1f}%)")

    # All clustering outputs go in a `clustering/` subdir next to the prefix, so
    # they don't clutter analysis/ (a run can produce dozens of rep PDBs). The
    # file basenames keep the prefix's name, e.g.
    #   analysis/remd_rep000  →  analysis/clustering/remd_rep000_cluster_*
    prefix_path = Path(args.prefix)
    out_dir = prefix_path.parent / "clustering"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_base = str(out_dir / prefix_path.name)

    write_assignments(f"{out_base}_cluster_assignments.csv",
                      frame_indices, times_ps, labels)
    reps = write_rep_pdbs(u, coords, labels, frame_indices, out_base)
    write_summary(f"{out_base}_cluster_summary.txt", args.method, n_frames,
                  labels, frame_indices, times_ps, reps,
                  args.cutoff, min_samples, args.n_clusters)
    plot_populations(f"{out_base}_cluster_populations.png", labels, n_frames)
    plot_timeseries(f"{out_base}_cluster_timeseries.png", times_ps, labels)
    print(f"[OK] clustering done → {out_base}_cluster_*")


if __name__ == "__main__":
    main()
