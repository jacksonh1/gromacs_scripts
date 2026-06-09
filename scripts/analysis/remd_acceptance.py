#!/usr/bin/env python3
"""
remd_acceptance.py — Exchange acceptance rates for a GROMACS T-REMD run.

Parses the "Replica exchange statistics" block from the GROMACS log and reports
per-pair empirical acceptance rates, mean Metropolis probabilities, and exchange
counts.

Usage:
    python remd_acceptance.py OUTDIR [--rep REP] [--plot]

    OUTDIR   path to the job output directory (contains prod/, analysis/)
    --rep    replica log to parse (default: 000; any log works — all replicas
             record the same Replica exchange statistics block)
    --plot   write a bar chart to OUTDIR/analysis/remd_acceptance.png
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def parse_header(log_text):
    m = re.search(r'There are (\d+) replicas', log_text)
    if not m:
        sys.exit("[ERROR] Could not find replica count in log.")
    n_replicas = int(m.group(1))

    m = re.search(r'Replica exchange interval:\s+(\d+)', log_text)
    if not m:
        sys.exit("[ERROR] Could not find exchange interval in log.")
    replex_interval = int(m.group(1))

    return {'n_replicas': n_replicas, 'replex_interval': replex_interval}


def parse_stats_block(log_text):
    m = re.search(r'(\d+) attempts, (\d+) odd, (\d+) even', log_text)
    if not m:
        sys.exit(
            "[ERROR] No 'Replica exchange statistics' block found.\n"
            "        Has the simulation completed?"
        )
    total = int(m.group(1))
    n_odd  = int(m.group(2))
    n_even = int(m.group(3))

    def parse_section(label):
        # Match: "Repl  <label>:\n" + index line + values line
        pat = re.compile(
            rf'Repl\s+{re.escape(label)}:\s*\nRepl[^\n]+\nRepl\s+([\d. ]+)',
            re.MULTILINE,
        )
        match = pat.search(log_text)
        if not match:
            sys.exit(f"[ERROR] Could not parse '{label}' section from log.")
        return [float(v) for v in match.group(1).split()]

    avg_prob = parse_section('average probabilities')
    n_exc    = parse_section('number of exchanges')
    avg_exc  = parse_section('average number of exchanges')

    # Pair i is attempted on even-type steps (n_even attempts) if i is even,
    # odd-type steps (n_odd attempts) if i is odd.
    attempts = [n_even if i % 2 == 0 else n_odd for i in range(len(n_exc))]

    return {
        'total_attempts': total,
        'n_odd':  n_odd,
        'n_even': n_even,
        'avg_prob':        avg_prob,
        'n_exchanges':     [int(x) for x in n_exc],
        'acceptance_rate': avg_exc,
        'attempts_per_pair': attempts,
    }


def get_temperatures(outdir, n_replicas):
    temps = []
    for i in range(n_replicas):
        log_path = Path(outdir) / 'prod' / f'rep{i:03d}' / 'remd.log'
        if not log_path.exists():
            sys.exit(f"[ERROR] Replica log not found: {log_path}")
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
            sys.exit(f"[ERROR] ensemble-temperature not found in {log_path}")
        temps.append(float(m.group(1)))
    return temps


def report(stats, temps, outdir, plot):
    n_pairs = len(stats['acceptance_rate'])
    rates   = stats['acceptance_rate']
    probs   = stats['avg_prob']
    n_exc   = stats['n_exchanges']
    atts    = stats['attempts_per_pair']

    # ── Console table ─────────────────────────────────────────────────────────
    hdr = (
        f"{'Pair':>8}  {'T_lo (K)':>9}  {'T_hi (K)':>9}"
        f"  {'Attempts':>9}  {'Exchanges':>9}  {'Rate':>8}  {'Avg. prob.':>10}"
    )
    print()
    print(
        f"REMD exchange acceptance rates"
        f"  ({stats['n_replicas']} replicas,"
        f" {stats['total_attempts']} attempts:"
        f" {stats['n_odd']} odd, {stats['n_even']} even)"
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

    log_path = Path(args.outdir) / 'prod' / f'rep{args.rep}' / 'remd.log'
    if not log_path.exists():
        sys.exit(f"[ERROR] Log not found: {log_path}")

    log_text = log_path.read_text(errors='replace')

    header = parse_header(log_text)
    stats  = parse_stats_block(log_text)
    stats['n_replicas']     = header['n_replicas']
    stats['replex_interval'] = header['replex_interval']

    temps = get_temperatures(args.outdir, header['n_replicas'])
    report(stats, temps, args.outdir, args.plot)


if __name__ == '__main__':
    main()
