# REST2 GPU exchange-energy noise — investigation and decision

**Status:** RESOLVED (2026-10-17). Measured, decomposed, and decided. Bottom line: on GPU,
REST2 (`-hrex`) exchange decisions carry a small, zero-mean energy noise. It does **not**
meaningfully bias the sampled ensemble for structural observables; **accept it** on the stock
GPU build. The `patch_gmx2023.5_gpu_energy_double.py` mitigation in the repo is **partial**
(fixes the minor source, not the dominant one) and is **not recommended** for production.

## The phenomenon

GROMACS evaluates the same configuration to a slightly different energy run-to-run on GPU:
the CUDA kernels reduce per-warp energies and `atomicAdd` them into single-precision (`float`)
accumulators, in a non-deterministic order. A `float` holding ~−3×10⁵ kJ/mol resolves only
~0.02 kJ/mol, so the running sum differs by ~1 kJ/mol between evaluations. CPU evaluation is
bit-exact.

REST2 replica exchange with PLUMED `-hrex` compares energies from **separate** evaluations
(PLUMED re-evaluates the Hamiltonian on the swapped coordinates, reported as `dplumed` in the
replica log; GROMACS zeroes its own `dE_term` under hrex by design). For two replicas with
**identical** Hamiltonians the exchange ΔE is analytically **zero**, so they must exchange with
probability 1 — but on GPU the two independent evaluations differ by the noise, so `dplumed` is
nonzero and acceptance falls below 100%. That shortfall is the diagnostic fingerprint.

Crucially, **forces are unaffected** — only reported energies. So each replica's continuous
dynamics at its own λ integrates the correct equations of motion and samples the correct
Boltzmann ensemble for that Hamiltonian. The noise only perturbs *which exchanges are accepted*.

## Measurements

All runs are identical-λ probes (two replicas, both λ = 1.0, `TEMPS_LIST` equivalent, 50
exchange attempts, `-replex 500`), on the `2023.5-plumed` build, L40S GPUs, unless noted.
`σ` is the standard deviation of `dplumed` (kT at 300 K); acceptance is the engine's own gate.

### Decomposition: two noise sources, PME dominates

System: `helix_fusion` built with **a99SB-disp** (4-site TIP4P-D water), ~15k atoms.

| build | PME | σ(dplumed) | acceptance | job |
|-------|-----|-----------|------------|-----|
| stock | `-pme gpu` | 0.42 kT | 0.81 | 25972453 |
| dblE (energy-double patch) | `-pme gpu` | 0.31 kT | 0.91 | 26161737 |
| stock | `-pme cpu` | 0.074 kT | 0.97 | 26201677 |
| dblE | `-pme cpu` | 0.0012 kT | 1.00 | 26161974 |

Reading the deltas:
- Moving **PME off the GPU** (stock, `-pme gpu` → `-pme cpu`): 0.42 → 0.074 kT. **PME
  reciprocal-space** GPU energy reduction is the **dominant** source.
- The **energy-double patch** (stock → dblE at `-pme gpu`): 0.42 → 0.31 kT. The short-range
  nonbonded (`nbnxm`) `float` accumulator — the only thing the patch touches — is the **minor**
  source.
- Both removed (dblE + `-pme cpu`): machine-exact (0.0012 kT, 100%).

A single-point determinism check (same frame re-evaluated 50×, `helix_fusion`, amber99sb-ildn,
stock GPU, job 25865807) independently confirmed the mechanism: **20 distinct total-potential
values in 50 GPU evaluations; 1 value in 20 on CPU** (bit-exact).

### On the real production system

Identical-λ probe on the actual TNFa-peptide REST2 system (**amber99sb-ildn / tip3p**, 20,959
atoms, stock GPU, job 26202725), i.e. exactly what production uses:

- **σ(dplumed) = 0.18 kT**, identical-λ acceptance 0.96.
- **Detailed-balance bias ≈ σ²/2 = 0.017 kT** (~1.6% of kT).

σ is smaller here than for the a99SB-disp helix because tip3p is **3-site**: fewer PME charge
sites → less PME reciprocal-space work → less PME noise. a99SB-disp's 4-site water is a noise
*amplifier* (and makes `-pme cpu` more expensive), but the phenomenon is force-field-agnostic.

### Production run is well-mixed despite the noise

The 100 ns production run (job 25087741: amber99sb-ildn/tip3p, 10 replicas, λ 1.0–0.75, 1 ps
exchange, NPT, stock GPU) shows the noise did not compromise sampling:

- acceptance **0.38–0.42** on every pair (healthy; un-depressed),
- **119–184 round trips per walker** over 100 ns (full-ladder traversal 100+×),
- slot-0 dwell ~0.03–0.21 (ideal 0.10 for 10 replicas; roughly uniform, no gross skew).

### Performance cost of the exact fix

`-pme cpu` is ~**3× slower** (helix, 2 replicas, 1 L40S: 369 → 122 ns/day). The energy-double
patch runs at full GPU speed (384 vs 369 ns/day) but, as above, only removes the minor source.

## Why the energy-double patch was partial, and why PME is hard to patch

The PR that added REST2's a99SB-disp support also shipped an energy-double GROMACS patch
(`patch_gmx2023.5_gpu_energy_double.py` + `gmx2023.5_gpu_energy_double.diff`). It made the
**nbnxm** energy accumulators (`eLJ`/`eElec`) double: correct, minimal (forces/trajectories
untouched), and it applied cleanly (10 edits, each matching once in the 2023.5 source). It was
built and validated here (job 26116676 built the `-dblE` variant; the identical-λ runs above
used it). **It was removed from the repo on 2026-10-17** (retrievable from git history) after
this investigation showed it addresses only the short-range source.

The dominant PME noise lives in `src/gromacs/ewald/pme_solve.cu`, where the reciprocal energy
is accumulated in a `float` per-thread sum and reduced via a packed warp-shuffle + shared-memory
reduction that bundles **7 components (6 virial + 1 energy)** into aligned float lanes, then
`atomicAdd`s into a global `float`. The same double-accumulator idea applies in principle, but
it is a far more invasive change (the float-specific packing must be rewritten) and risks
corrupting the **virial** (hence pressure coupling). It was not attempted.

## Decision

**Accept the stock-GPU bias for the lab's REST2 use (structural characterization of designed
binders: RMSD, RMSF, secondary structure, conformational populations).** Justification:

- Forces are exact → each λ samples the correct Boltzmann ensemble.
- The noise is zero-mean → no systematic energy offset.
- The detailed-balance bias is second-order and tiny: **~0.017 kT** (amber99sb-ildn/tip3p),
  ~0.09 kT for a99SB-disp — far below statistical error and the force field's own ~kT accuracy.
- The real production run mixed thoroughly (100+ round trips).

Use **`-pme cpu`** (exact; ~3× slower) only if a future use case needs high-precision free
energies or runs a very fine ladder with marginal adjacent-pair acceptance.

**The energy-double patch is not the answer** (and was removed): it is a self-maintained GROMACS
fork that removes only the minor noise source, leaving ~0.31 kT / 91% at full GPU speed — no
practical benefit over the stock build. If a future need justifies a fast-and-exact GPU path,
the real work is doubling the PME solve-kernel energy reduction (invasive, virial-risky), not
the nbnxm accumulator.

## Safeguards (already in the pipeline)

- **Identical-λ / scale-1.0 acceptance self-check** — the real failure mode to guard against is
  a *broken* hrex build (the 2024.3 build gives zero exchanges silently). `check_hrex_acceptance.py`
  catches that; the small GPU noise is not that failure.
- **Round-trip / convergence check** (`gromd-roundtrip`, observable stability) — confirms enough
  sampling was achieved despite the slightly reduced acceptance.

## Reproduction

Identical-λ probe: build (or reuse a run's) solvated system, make the λ = 1.0 topology
(`mark_hot_region.py` + `plumed partial_tempering 1.0`), grompp two replicas from it with
different `gen-seed`, and run the engine's exact command
(`mpirun -np 2 gmx_mpi mdrun -multidir rep000 rep001 -deffnm rest2 -plumed plumed.dat -hrex
-replex 500 -nsteps 25000 -nb gpu -pme gpu -bonded gpu`). Parse `dplumed` from `rep000/rest2.log`;
σ is the per-exchange noise, σ²/2 the detailed-balance bias. Swap `-pme gpu`→`-pme cpu` or the
stock build for the `-dblE` build to isolate the two sources.
