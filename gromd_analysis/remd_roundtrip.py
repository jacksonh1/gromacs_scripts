#!/usr/bin/env python3
"""
remd_roundtrip.py — Replica mixing diagnostics for a GROMACS T-REMD / REST2 run.

Acceptance rate (gromd-acceptance) is a *local* quantity: the probability a
neighbouring pair swaps. It does not tell you whether configurations actually
traverse the whole temperature ladder. This tool answers the *global* mixing
question by reconstructing the permutation walk from the per-attempt `Repl ex`
lines in the log and reporting, per configuration ("walker"):

  * round trips        — full T_min -> T_max -> T_min cycles completed
  * round-trip time     — run length / mean round trips per walker
  * temperature dwell  — fraction of the run each walker spent in each slot

Why this matters: the constant-temperature ensemble at a fixed slot (e.g.
rep000 = T_min) is *unbiased* regardless of mixing — every frame there is a
valid Boltzmann sample by construction. What poor mixing costs is the *number
of independent samples* that slot receives: each round trip carries a walker up
to T_max (where it crosses barriers / decorrelates) and back down, delivering
one effectively fresh configuration to T_min. Few round trips => large
statistical error bars, not a wrong mean. Round-trip time is therefore the
efficiency metric; acceptance alone can look healthy while the ladder still
fails to mix if a single pair bottlenecks.

Correctness gate: the swap parsing is validated against GROMACS's own numbers.
We recount, from the `Repl ex` lines we parse, how many times each neighbouring
pair exchanged, and assert it equals the `number of exchanges` list in the
run's `Replica exchange statistics` block (parsed by the shared, tested
ExchangeStats). A mismatch means our permutation reconstruction disagrees with
GROMACS and every downstream number would be wrong, so we raise instead.

Round-trip definition (bottom-anchored, the standard one): a walker completes a
round trip each time it arrives at the bottom slot (T_min, slot 0) having
visited the top slot (T_max, slot N-1) at least once since its previous arrival
at the bottom. This counts full down-and-back cycles and is what the round-trip
*time* (run length / round trips) is defined against. It deliberately does not
credit partial excursions (bottom -> middle -> bottom) or a final leg that
reaches the top but never returns.

Usage:
    gromd-roundtrip OUTDIR [--rep REP] [--plot]

    OUTDIR   path to the job output directory (contains prod/, analysis/)
    --rep    replica log to parse (default: 000; every replica records the same
             global `Repl ex` decisions, so any one log reconstructs the walk)
    --plot   write bar charts to OUTDIR/analysis/remd_roundtrips.png
"""

import argparse
import csv
import re
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path

from gromd_analysis.remd_log import (
    ExchangeStats,
    RemdLogError,
    get_temperatures,
    prod_basename,
)


@dataclass(frozen=True)
class MixingStats:
    """Reconstructed replica-mixing diagnostics.

    All per-walker lists have length `n_replicas` and are indexed by walker id
    (the slot the walker started in). `dwell_frac[w]` is a length-`n_replicas`
    list: the fraction of exchange attempts walker `w` spent in each slot.
    """

    n_replicas: int
    n_attempts: int               # number of Repl ex records parsed
    run_ns: float                 # wall time of the last exchange, in ns
    round_trips: list[int]        # per walker: completed T_min<->T_max cycles
    dwell_frac: list[list[float]] # per walker: occupancy fraction per slot

    @property
    def mean_round_trips(self) -> float:
        return sum(self.round_trips) / self.n_replicas

    @property
    def round_trip_time_ns(self) -> float:
        """Run length / mean round trips per walker. inf if nothing round-tripped."""
        m = self.mean_round_trips
        return self.run_ns / m if m > 0 else float('inf')


def _parse_replex_line(line: str, n_replicas: int) -> list[int]:
    """Return the list of lower-slot indices that swapped on this attempt.

    A `Repl ex` line prints the slot axis 0,1,...,N-1 in order with an `x`
    between two adjacent slot columns that exchanged, e.g.

        Repl ex  0    1    2    3 x  4    5 x  6 ...

    means slots (3,4) and (5,6) swapped. We walk the tokens after `Repl ex`,
    counting numeric columns to know the current slot, and record `cur` for each
    `x` (the swap is between `cur` and `cur+1`). Swaps within one attempt are
    always non-overlapping (GROMACS alternates odd/even pairs), so the order of
    application does not matter.
    """
    toks = line.split()[2:]  # drop the "Repl" "ex" prefix
    cur = -1
    swaps = []
    for t in toks:
        if t == 'x':
            swaps.append(cur)
        else:
            cur += 1
    # The line must enumerate exactly n_replicas slot columns; a short/long line
    # would silently shift every swap index, so fail loudly instead.
    if cur != n_replicas - 1:
        raise RemdLogError(
            f"'Repl ex' line lists {cur + 1} slot columns, expected {n_replicas}: {line!r}")
    return swaps


def _read_log(log_path: Path):
    """One streaming pass over the (possibly multi-GB) log.

    Returns (replex_lines, last_time_ps, stats_text). `stats_text` is the head
    plus tail of the log — enough for the shared ExchangeStats parser (replica
    count + interval live in the setup header; the statistics block is at the
    very end) without holding the whole file in memory.
    """
    replex_lines = []
    last_time_ps = 0.0
    head = []
    tail = deque(maxlen=400)
    with open(log_path, errors='replace') as f:
        for lineno, line in enumerate(f):
            if line.startswith("Repl ex"):
                replex_lines.append(line)
            elif line.startswith("Replica exchange at step"):
                m = re.search(r'time\s+([\d.]+)', line)
                if m:
                    last_time_ps = float(m.group(1))
            if lineno < 600:
                head.append(line)
            tail.append(line)
    if not replex_lines:
        raise RemdLogError(
            f"no 'Repl ex' lines in {log_path}; not a replica-exchange run, or unfinished.")
    return replex_lines, last_time_ps, ''.join(head) + ''.join(tail)


def compute(log_path: Path) -> MixingStats:
    replex_lines, last_time_ps, stats_text = _read_log(log_path)

    # Ground truth for the self-check, via the shared, tested parser.
    stats = ExchangeStats.parse(stats_text)
    n = stats.n_replicas

    # Reconstruct the permutation walk. occ[slot] = walker currently in that slot.
    occ = list(range(n))
    pair_swaps = [0] * (n - 1)
    dwell = [[0] * n for _ in range(n)]  # dwell[walker][slot]
    top = n - 1
    seen_top = [False] * n
    round_trips = [0] * n

    for line in replex_lines:
        for lo in _parse_replex_line(line, n):
            pair_swaps[lo] += 1
            occ[lo], occ[lo + 1] = occ[lo + 1], occ[lo]
        # occ must stay a permutation; a duplicate would corrupt every count.
        assert len(set(occ)) == n, "reconstructed occupancy is not a permutation"
        for slot, w in enumerate(occ):
            dwell[w][slot] += 1
            if slot == top:
                seen_top[w] = True
            elif slot == 0 and seen_top[w]:
                round_trips[w] += 1
                seen_top[w] = False

    # ── CORRECTNESS GATE ─────────────────────────────────────────────────────
    # Our parsed swaps, recounted per pair, must equal GROMACS's own tally.
    if pair_swaps != stats.n_exchanges:
        diffs = [(i, pair_swaps[i], stats.n_exchanges[i])
                 for i in range(n - 1) if pair_swaps[i] != stats.n_exchanges[i]]
        raise RemdLogError(
            "self-check FAILED: reconstructed per-pair exchange counts disagree with "
            "the log's 'number of exchanges' block — the permutation walk is wrong, "
            f"so round trips would be too. Mismatches (pair, recon, gmx): {diffs}")

    n_attempts = len(replex_lines)
    dwell_frac = [[c / n_attempts for c in row] for row in dwell]

    return MixingStats(
        n_replicas=n,
        n_attempts=n_attempts,
        run_ns=last_time_ps / 1000.0,
        round_trips=round_trips,
        dwell_frac=dwell_frac,
    )


def report(mix: MixingStats, temps, outdir, plot):
    n = mix.n_replicas
    top = n - 1
    ideal = 1.0 / n

    print()
    print(
        f"REMD replica mixing"
        f"  ({n} replicas, {mix.n_attempts} exchange attempts, {mix.run_ns:.1f} ns)"
    )
    # For T-REMD these differ (a real temperature gradient); for REST2 every
    # replica reports 300 K (the ladder is in effective temperature), so the
    # endpoints are named by slot, with the reported temperature in parentheses.
    print(f"  endpoints: slot 0 ({temps[0]:.1f} K)  <->  slot {top} ({temps[top]:.1f} K)")
    print(f"  [self-check PASSED: per-pair swaps match the log's exchange counts]")
    print()

    hdr = f"{'Walker':>7}  {'Round trips':>11}  {'RT time (ns)':>12}  {'slot0 dwell':>11}"
    print(hdr)
    print("-" * len(hdr))
    for w in range(n):
        rt = mix.round_trips[w]
        rt_time = f"{mix.run_ns / rt:.1f}" if rt > 0 else "inf"
        print(
            f"{w:>7}  {rt:>11d}  {rt_time:>12}  {mix.dwell_frac[w][0]:>11.3f}")

    print()
    print(
        f"Mean round trips/walker: {mix.mean_round_trips:.1f}"
        f"   Min: {min(mix.round_trips)}"
        f"   Round-trip time: {mix.round_trip_time_ns:.2f} ns"
    )
    print(
        f"slot0 occupancy spread: {min(r[0] for r in mix.dwell_frac):.3f}"
        f"–{max(r[0] for r in mix.dwell_frac):.3f}   (ideal uniform = {ideal:.3f})"
    )
    print(
        "Note: round-trip time << run length means the ladder mixes well; each round\n"
        "      trip delivers one decorrelated configuration to slot 0. Dwell imbalance\n"
        "      over a finite run is scatter, not ensemble bias — the slot-0 ensemble\n"
        "      stays a correct Boltzmann ensemble regardless."
    )
    print()

    analysis_dir = Path(outdir) / 'analysis'
    analysis_dir.mkdir(exist_ok=True)
    csv_path = analysis_dir / 'remd_roundtrips.csv'
    with open(csv_path, 'w', newline='') as f:
        wtr = csv.writer(f)
        wtr.writerow(['walker', 'round_trips', 'roundtrip_time_ns', 'slot0_dwell_frac'])
        for w in range(n):
            rt = mix.round_trips[w]
            rt_time = f'{mix.run_ns / rt:.3f}' if rt > 0 else 'inf'
            wtr.writerow([w, rt, rt_time, f'{mix.dwell_frac[w][0]:.4f}'])
    print(f"[OK] CSV written to: {csv_path}")

    if plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("[WARN] matplotlib not available; skipping plot.")
            return
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(9, n * 0.7), 4))
        ax1.bar(range(n), mix.round_trips, color='steelblue', alpha=0.8)
        ax1.axhline(mix.mean_round_trips, color='crimson', linestyle='--',
                    linewidth=0.9, label=f'mean {mix.mean_round_trips:.1f}')
        ax1.set_xlabel('Walker (starting slot)')
        ax1.set_ylabel('Round trips')
        ax1.set_title(f'Round trips  (τ = {mix.round_trip_time_ns:.1f} ns)')
        ax1.legend()
        ax2.bar(range(n), [r[0] * 100 for r in mix.dwell_frac],
                color='seagreen', alpha=0.8)
        ax2.axhline(ideal * 100, color='crimson', linestyle='--',
                    linewidth=0.9, label=f'uniform {ideal*100:.1f}%')
        ax2.set_xlabel('Walker (starting slot)')
        ax2.set_ylabel(f'Time in slot 0 ({temps[0]:.0f} K) (%)')
        ax2.set_title('slot-0 occupancy per walker')
        ax2.legend()
        fig.tight_layout()
        png_path = analysis_dir / 'remd_roundtrips.png'
        fig.savefig(png_path, dpi=150)
        print(f"[OK] Plot written to:  {png_path}")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description='Report REMD/REST2 replica mixing (round trips, dwell) from a GROMACS log.'
    )
    ap.add_argument('outdir', help='job output directory')
    ap.add_argument('--rep', default='000',
                    help='replica log to parse (default: 000; any replica log works)')
    ap.add_argument('--plot', action='store_true',
                    help='write mixing bar charts to OUTDIR/analysis/remd_roundtrips.png')
    args = ap.parse_args()

    try:
        basename = prod_basename(args.outdir)
        log_path = Path(args.outdir) / 'prod' / f'rep{args.rep}' / f'{basename}.log'
        if not log_path.exists():
            raise RemdLogError(f"Log not found: {log_path}")

        mix = compute(log_path)
        temps = get_temperatures(args.outdir, mix.n_replicas, basename)
        report(mix, temps, args.outdir, args.plot)
    except RemdLogError as exc:
        sys.exit(f"[ERROR] {exc}")


if __name__ == '__main__':
    main()
