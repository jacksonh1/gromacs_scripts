# GROMACS T-REMD Output Reference

Notes on the structure and content of GROMACS T-REMD output files, accumulated
while building the analysis pipeline. Use this before parsing any log or trajectory.

---

## File layout

```
OUTDIR/
  prod/
    rep000/           ← temperature slot 0 (lowest T, e.g. 300 K)
      remd.tpr        ← run input file (topology + positions + velocities)
      remd.log        ← main log: mdp params, per-exchange records, end-of-run stats
      remd.xtc        ← trajectory at this temperature slot
      remd.edr        ← energy file at this temperature slot
    rep001/           ← temperature slot 1
    ...
    rep047/           ← temperature slot N-1 (highest T)
  trajectories/
    remd_rep000.xtc   ← concatenated/copied trajectory for convenient access
    ...
  analysis/           ← created by analysis scripts
```

**Critical:** each `rep{i}/` directory corresponds to a **fixed temperature slot**,
not a specific configuration (molecule). Coordinates are exchanged between slots at
each exchange step, so the molecule in `rep000/remd.xtc` at time t is not necessarily
the same physical configuration as at time t+dt.

---

## Temperature of each slot

Each replica log records its temperature in the mdp parameters section near the top:

```
ensemble-temperature           = 300
```

Read this line from `prod/rep{i:03d}/remd.log` to get the temperature for slot i.
The temperature section appears within the first ~200 lines of the log.

---

## Per-exchange blocks

Appears repeatedly throughout the log, one block per exchange attempt:

```
Replica exchange at step 500 time 1.00000
Repl 10 <-> 11  dE_term =  2.208e-01 (kT)
Repl ex  0    1    2    3    4 x  5    6 x  7    8 x  9   10 x 11 ...
Repl pr   .08       .04       1.0       .17       .64 ...
```

- `step N` — simulation step at which this exchange was attempted
- `Repl X <-> Y  dE_term` — the pair that was evaluated for exchange (only one pair
  printed here, but multiple pairs are attempted; this line reflects the "worst case"
  or the pair that the log file's replica is directly involved in)
- `Repl ex` — full ladder showing which pairs exchanged: `x` between two adjacent
  replica indices means that pair accepted the exchange
- `Repl pr` — mean Metropolis acceptance probabilities for each attempted pair,
  in the same order as the attempted pairs

### Alternating even/odd pattern

GROMACS alternates which pairs are attempted each exchange step:

```
exchange_number = step // replex_interval      # (1, 2, 3, ...)
if exchange_number is odd:   attempt even pairs  (0–1, 2–3, 4–5, ...)
if exchange_number is even:  attempt odd pairs   (1–2, 3–4, 5–6, ...)
```

Confirmed: step 500 (exchange_number=1, odd) → even pairs; step 1000 (exchange_number=2, even) → odd pairs.

For each pair, it is attempted on either `n_even` or `n_odd` steps (see statistics
block below), depending on whether it is an even- or odd-indexed pair.

### Parsing `Repl ex` for accepted pairs

```python
tokens = line[len("Repl ex"):].split()
for k, tok in enumerate(tokens):
    if tok == 'x':
        left, right = int(tokens[k-1]), int(tokens[k+1])
        # pair (left, right) exchanged at this step
```

---

## End-of-run statistics block

GROMACS pre-computes acceptance statistics and writes them once at the end of the log.
**Use this block — do not reparse the per-exchange lines.**

```
Replica exchange statistics
Repl  1999 attempts, 1000 odd, 999 even
Repl  average probabilities:
Repl     0    1    2    3  ...  46   47
Repl      .48  .47  .51  .52 ...  .57  .55
Repl  number of exchanges:
Repl     0    1    2    3  ...  46   47
Repl      474  457  510  517 ...  561  561
Repl  average number of exchanges:
Repl     0    1    2    3  ...  46   47
Repl      .47  .46  .51  .52 ...  .56  .56
```

- **Column index i** = pair i ↔ i+1 (0-indexed); N replicas → N-1 values
- **`average probabilities`** = mean Metropolis criterion: E[min(1, exp(−ΔE/kT))]
  (theoretical average, not empirical rate)
- **`number of exchanges`** = count of accepted exchanges per pair
- **`average number of exchanges`** = empirical acceptance rate
  (= exchanges / attempts_for_that_pair, already accounting for even/odd alternation)
- **Attempts per pair**: even-indexed pair i uses `n_even` attempts; odd-indexed uses `n_odd`

All replica logs record the same statistics block (synchronized via MPI); any one log
can be used to extract it.

### Parsing the statistics block

```python
import re

m = re.search(r'(\d+) attempts, (\d+) odd, (\d+) even', log_text)
total, n_odd, n_even = int(m.group(1)), int(m.group(2)), int(m.group(3))

def parse_section(label):
    pat = re.compile(
        rf'Repl\s+{re.escape(label)}:\s*\nRepl[^\n]+\nRepl\s+([\d. ]+)',
        re.MULTILINE,
    )
    return [float(v) for v in pat.search(log_text).group(1).split()]

avg_prob = parse_section('average probabilities')
n_exc    = parse_section('number of exchanges')
avg_exc  = parse_section('average number of exchanges')
```

---

## Empirical Transition Matrix

Follows the statistics block:

```
Repl       1       2  ...     48
Repl  0.7629  0.2371  0.0000 ...  0
Repl  0.2371  0.5343  0.2286 ...  1
```

- Column headers are 1-indexed (1..N); row ends with 0-indexed configuration index
- Entry (i, j) = empirical probability of transitioning from slot i to slot j in one
  exchange step
- **This is NOT dwell time.** It is the one-step transition probability matrix.
  The diagonal element P(i→i) = 1 − (acceptance rate for pairs adjacent to i × 0.5),
  approximately.

---

## Dwell time and round trips

**These are not directly available from any single output file.**

Computing them requires reconstructing which physical configuration is at which
temperature slot at each exchange step — i.e., tracking a permutation array that
is updated with each `Repl ex` line. Since each log corresponds to a temperature
slot (not a configuration), this reconstruction requires parsing all per-exchange
`Repl ex` lines from one log and applying them cumulatively.

This is feasible but non-trivial; deferred to a future script.

---

## Exchange acceptance rate targets

| Rate      | Interpretation                                              |
|-----------|-------------------------------------------------------------|
| < 15%     | Replicas too far apart; poor mixing                        |
| 20–30%    | Well-tuned T-REMD (typical recommendation)                 |
| > 50%     | Replicas too close together; wasting compute resources      |

A uniform acceptance rate across the ladder means temperature spacing is even.
Rising rates toward higher T (as seen in the example job) suggest the spacing
could be compressed at low T and expanded at high T.

---

## PBC correction for REMD trajectories

Use `-pbc mol`, NOT `-pbc nojump`. The `nojump` algorithm compares consecutive
frames to detect box-boundary crossings, but T-REMD coordinate exchanges cause
discontinuous jumps between frames that `nojump` misinterprets. `-pbc mol` is a
per-frame operation and is safe.

```bash
printf "Protein\nSystem\n" | gmx trjconv \
  -s prod/rep000/remd.tpr \
  -f trajectories/remd_rep000.xtc \
  -o analysis/remd_rep000_pbc.xtc \
  -pbc mol -center -ur compact
```

---

## `gmx trjconv` atom-count matching

The `-s` reference passed to `gmx trjconv` must have the **same atom count** as
the XTC. Passing a protein-only GRO as `-s` while the XTC contains the full system
causes GROMACS to silently truncate the trajectory read, producing a degenerate
fitting matrix and a fatal "Too many iterations in routine JACOBI" error.

When stripping solvent and aligning in one call, extract a full-system GRO as the
`-s` reference and a protein-only GRO as the user-facing reference separately:

```bash
printf "System\n"  | gmx trjconv -s remd.tpr -f pbc.xtc -o ref_full.gro -dump 0
printf "Protein\n" | gmx trjconv -s remd.tpr -f pbc.xtc -o ref.gro      -dump 0
printf "Backbone\nProtein\n" | gmx trjconv \
  -s ref_full.gro -f pbc.xtc -o fit.xtc -fit rot+trans
```
