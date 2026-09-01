#!/usr/bin/env python3
"""
remd_log.py — Exchange acceptance rates for a GROMACS T-REMD run.

Parses the "Replica exchange statistics" block from the GROMACS log and reports
per-pair empirical acceptance rates, mean Metropolis probabilities, and exchange
counts.

Usage:
    gromd-acceptance OUTDIR [--rep REP] [--plot]

    OUTDIR   path to the job output directory (contains prod/, analysis/)
    --rep    replica log to parse (default: 000; any log works — all replicas
             record the same Replica exchange statistics block)
    --plot   write a bar chart to OUTDIR/analysis/remd_acceptance.png
"""

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


class RemdLogError(Exception):
    """The log is not a finished replica-exchange run, or its statistics block is unparseable."""


@dataclass(frozen=True)
class ExchangeStats:
    """The `Replica exchange statistics` block, parsed.

    GROMACS pre-computes these at the end of the run — per-pair acceptance rates
    and counts are read straight out of the block, never recounted from per-frame
    `Repl ex` lines.

    The four per-pair lists are parallel and all have length `n_replicas - 1`
    (one entry per neighbouring temperature pair): pair `i` is slot `i` ↔ `i+1`.
    """

    n_replicas: int
    replex_interval: int          # exchange attempt interval, in steps
    total_attempts: int
    n_odd: int                    # attempts on odd-type steps
    n_even: int                   # attempts on even-type steps
    avg_prob: list[float]         # per pair: mean Metropolis probability
    n_exchanges: list[int]        # per pair: accepted exchanges
    acceptance_rate: list[float]  # per pair: empirical rate (avg number of exchanges)
    attempts_per_pair: list[int]  # per pair: attempts of that pair's parity

    @property
    def n_pairs(self) -> int:
        return len(self.acceptance_rate)

    @classmethod
    def parse(cls, log_text: str) -> "ExchangeStats":
        """Parse a replica log, or raise RemdLogError. Any replica's log will do —
        every replica records the same statistics block."""
        n_replicas = _search_int(
            r'There are (\d+) replicas', log_text, "replica count")
        replex_interval = _search_int(
            r'Replica exchange interval:\s+(\d+)', log_text, "exchange interval")

        m = re.search(r'(\d+) attempts, (\d+) odd, (\d+) even', log_text)
        if not m:
            raise RemdLogError(
                "No 'Replica exchange statistics' block found. Has the simulation completed?")
        total, n_odd, n_even = (int(g) for g in m.groups())

        avg_prob = _parse_section('average probabilities', log_text)
        n_exc    = _parse_section('number of exchanges', log_text)
        avg_exc  = _parse_section('average number of exchanges', log_text)

        # N replicas ⇒ N-1 neighbouring pairs. A mismatch means the block was
        # truncated or the format changed; either way the per-pair table below
        # would be silently wrong, so fail here instead.
        for label, values in (('average probabilities', avg_prob),
                              ('number of exchanges', n_exc),
                              ('average number of exchanges', avg_exc)):
            if len(values) != n_replicas - 1:
                raise RemdLogError(
                    f"'{label}' has {len(values)} values, expected {n_replicas - 1} "
                    f"for {n_replicas} replicas.")

        # Pair i is attempted on even-type steps (n_even attempts) if i is even,
        # odd-type steps (n_odd attempts) if i is odd.
        attempts = [n_even if i % 2 == 0 else n_odd for i in range(len(n_exc))]

        return cls(
            n_replicas=n_replicas,
            replex_interval=replex_interval,
            total_attempts=total,
            n_odd=n_odd,
            n_even=n_even,
            avg_prob=avg_prob,
            n_exchanges=[int(x) for x in n_exc],
            acceptance_rate=avg_exc,
            attempts_per_pair=attempts,
        )


def _search_int(pattern: str, log_text: str, what: str) -> int:
    m = re.search(pattern, log_text)
    if not m:
        raise RemdLogError(f"Could not find {what} in log.")
    return int(m.group(1))


def _parse_section(label: str, log_text: str) -> list[float]:
    """Values from one `Repl  <label>:` section: label line, index line, value line."""
    pat = re.compile(
        rf'Repl\s+{re.escape(label)}:\s*\nRepl[^\n]+\nRepl\s+([\d. ]+)',
        re.MULTILINE,
    )
    match = pat.search(log_text)
    if not match:
        raise RemdLogError(f"Could not parse '{label}' section from log.")
    return [float(v) for v in match.group(1).split()]


def prod_basename(outdir):
    """Production log basename: 'remd' for T-REMD, 'rest2' for REST2 (same log format
    and exchange-statistics block; REST2 replicas all report the same physical
    temperature, so the per-pair acceptance rates are the meaningful output)."""
    for base in ('remd', 'rest2'):
        if (Path(outdir) / 'prod' / 'rep000' / f'{base}.log').exists():
            return base
    raise RemdLogError(f"no prod/rep000/{{remd,rest2}}.log under {outdir}")


def get_temperatures(outdir, n_replicas, basename):
    temps = []
    for i in range(n_replicas):
        log_path = Path(outdir) / 'prod' / f'rep{i:03d}' / f'{basename}.log'
        if not log_path.exists():
            raise RemdLogError(f"Replica log not found: {log_path}")
        # Temperature appears in the mdp parameters section near the top of the log.
        # Read in chunks to avoid loading multi-GB files fully into memory.
        header_text = []
        with open(log_path, 'r', errors='replace') as f:
            for lineno, line in enumerate(f):
                header_text.append(line)
                if lineno > 300:
                    break
        text = ''.join(header_text)
        m = re.search(r'ensemble-temperature\s*=\s*([\d.]+)', text)
        if not m:
            raise RemdLogError(f"ensemble-temperature not found in {log_path}")
        temps.append(float(m.group(1)))
    return temps


def report(stats: ExchangeStats, temps, outdir, plot):
    # ASSUMES: temps is per-slot, so the pair table can index temps[i] and temps[i+1].
    if len(temps) != stats.n_replicas:
        raise RemdLogError(
            f"got {len(temps)} temperatures for {stats.n_replicas} replicas.")

    n_pairs = stats.n_pairs
    rates   = stats.acceptance_rate
    probs   = stats.avg_prob
    n_exc   = stats.n_exchanges
    atts    = stats.attempts_per_pair

    # ── Console table ─────────────────────────────────────────────────────────
    hdr = (
        f"{'Pair':>8}  {'T_lo (K)':>9}  {'T_hi (K)':>9}"
        f"  {'Attempts':>9}  {'Exchanges':>9}  {'Rate':>8}  {'Avg. prob.':>10}"
    )
    print()
    print(
        f"REMD exchange acceptance rates"
        f"  ({stats.n_replicas} replicas,"
        f" {stats.total_attempts} attempts:"
        f" {stats.n_odd} odd, {stats.n_even} even)"
    )
    print(hdr)
    print("-" * len(hdr))
    for i in range(n_pairs):
        pair_label = f"{i}–{i+1}"
        print(
            f"{pair_label:>8}  {temps[i]:>9.1f}  {temps[i+1]:>9.1f}"
            f"  {atts[i]:>9d}  {n_exc[i]:>9d}"
            f"  {rates[i]*100:>7.1f}%  {probs[i]*100:>9.0f}%"
        )

    mean_rate = sum(rates) / len(rates)
    min_rate  = min(rates)
    max_rate  = max(rates)
    print()
    print(
        f"Mean acceptance rate: {mean_rate*100:.1f}%"
        f"   Min: {min_rate*100:.1f}%"
        f"   Max: {max_rate*100:.1f}%"
    )
    print("Note: 20–30% per pair is typical for well-tuned T-REMD.")
    print()

    # ── CSV ───────────────────────────────────────────────────────────────────
    analysis_dir = Path(outdir) / 'analysis'
    analysis_dir.mkdir(exist_ok=True)
    csv_path = analysis_dir / 'remd_acceptance.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([
            'pair_lo', 'pair_hi', 'T_lo_K', 'T_hi_K',
            'attempts', 'exchanges', 'acceptance_rate', 'avg_metropolis_prob',
        ])
        for i in range(n_pairs):
            w.writerow([
                i, i + 1,
                f'{temps[i]:.3f}', f'{temps[i+1]:.3f}',
                atts[i], n_exc[i],
                f'{rates[i]:.4f}', f'{probs[i]:.4f}',
            ])
    print(f"[OK] CSV written to: {csv_path}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    if plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("[WARN] matplotlib not available; skipping plot.")
            return

        pair_labels = [f"{temps[i]:.0f}–{temps[i+1]:.0f}" for i in range(n_pairs)]
        fig, ax = plt.subplots(figsize=(max(8, n_pairs * 0.4), 4))
        ax.bar(range(n_pairs), [r * 100 for r in rates], color='steelblue', alpha=0.8)
        ax.axhline(20, color='green',  linestyle='--', linewidth=0.9, label='20%')
        ax.axhline(30, color='orange', linestyle='--', linewidth=0.9, label='30%')
        ax.set_xticks(range(n_pairs))
        ax.set_xticklabels(pair_labels, rotation=90, fontsize=6)
        ax.set_xlabel('Temperature pair (K)')
        ax.set_ylabel('Acceptance rate (%)')
        ax.set_title('REMD exchange acceptance rates')
        ax.legend(title='Target range')
        fig.tight_layout()
        png_path = analysis_dir / 'remd_acceptance.png'
        fig.savefig(png_path, dpi=150)
        print(f"[OK] Plot written to:  {png_path}")
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description='Report REMD exchange acceptance rates from a GROMACS log.'
    )
    ap.add_argument('outdir', help='job output directory')
    ap.add_argument(
        '--rep', default='000',
        help='replica log to parse (default: 000; any replica log works)',
    )
    ap.add_argument(
        '--plot', action='store_true',
        help='write acceptance rate bar chart to OUTDIR/analysis/remd_acceptance.png',
    )
    args = ap.parse_args()

    try:
        basename = prod_basename(args.outdir)
        log_path = Path(args.outdir) / 'prod' / f'rep{args.rep}' / f'{basename}.log'
        if not log_path.exists():
            raise RemdLogError(f"Log not found: {log_path}")

        stats = ExchangeStats.parse(log_path.read_text(errors='replace'))
        temps = get_temperatures(args.outdir, stats.n_replicas, basename)
        report(stats, temps, args.outdir, args.plot)
    except RemdLogError as exc:
        sys.exit(f"[ERROR] {exc}")


if __name__ == '__main__':
    main()
