# Task: validate a REST2 setup before production

You are validating a replica-exchange-with-solute-tempering (REST2) setup in OpenMM.
Work through the gates in order. Each gate states what to run, the numeric pass
criterion, and what it catches. **Stop at the first failed BLOCKING gate and report** —
do not continue to later gates and do not start production.

Report results only as the table at the bottom. Do not summarise a gate as "looks
reasonable": every gate has a number, so report the number.

## What REST2 is, so you can check it

Replica *m* runs at a single temperature `T_0` on the Hamiltonian

```
E_m = λ_m·E_pp + √λ_m·E_pw + E_ww        λ_m = T_0/T_m^eff,   λ ∈ (0, 1]
```

`E_pp` = intra-solute, `E_pw` = solute–solvent, `E_ww` = solvent–solvent. Implemented by
scaling, **for solute particles only**:

| Quantity | Factor | Why |
| :- | :- | :- |
| Charges | `√λ` | makes solute–solute electrostatics `λ` and solute–solvent `√λ`, exactly, PME included (the Ewald energy is bilinear in charge) |
| LJ `epsilon` | `λ` | same split via the geometric-mean combination rule |
| Solute torsions | `λ` | `E_pp` is scaled as a whole |
| CMAP grids | `λ` | a pure function of solute backbone dihedrals, so it is in `E_pp` |
| Pair-specific / off-diagonal LJ | `λ` solute–solute, `√λ` solute–solvent | not reachable by scaling per-particle `epsilon` |
| `NonbondedForce` exceptions (1-4) | `chargeProd`·λ, `epsilon`·λ when both particles are solute | stored explicitly, easy to miss |
| Bonds, angles | `1` | the REST2 paper reports scaling them does not help sampling |
| Everything involving only solvent | `1` | — |

Reference: Wang, Friesner & Berne, *J. Phys. Chem. B* 115, 9431 (2011),
doi:10.1021/jp204407d.

## Forcefield-specific traps to check for explicitly

- **CHARMM36m**: the CMAP grid is the trap. Scale every map:
  `CMAPTorsionForce.getMapParameters(i)` → multiply the energy vector by λ →
  `setMapParameters(i, size, scaled)`. If the solute is the only CMAP-bearing molecule,
  scaling all maps is equivalent to scaling the solute's.
- **a99SB-disp**: the trap is a *pair-specific* LJ modification between backbone carbonyl
  oxygens and backbone amide hydrogens, added to strengthen backbone hydrogen bonds
  (doi:10.1073/pnas.1800690115). Find how your forcefield XML implements it — a dedicated
  type pair, an exception, or a `CustomNonbondedForce` — by inspecting the built `System`,
  not the XML. Both atoms are solute, so the factor is λ.

In both forcefields the trap is a backbone secondary-structure driver. Left unscaled, hot
replicas keep full-strength backbone structure while sterics and sidechains melt.

## Setup required before gate 1

Put each force in its own force group (`force.setForceGroup(i)`) and record the mapping.
Build a **separate `System` per λ** rather than mutating one in place. Run all
single-point energy gates in double precision (`Platform.getPlatformByName('Reference')`,
or CUDA with `precision='double'`); the tolerances below are meaningless in mixed
precision. Use ≥20 frames spanning distinct conformations, not one frame.

---

## GATE 0 — System audit · BLOCKING · seconds

1. Enumerate `system.getForces()` for the λ=1 system. Classify every force as scaled-by-λ,
   invariant, or **unhandled**.
2. Across all λ systems, compare: particle count, every particle mass, constraint count
   and every constraint distance, bond and angle parameters.
3. Write the solute (hot) particle index set to a file and report its size.

**Pass**: zero unhandled forces; zero differences in masses or constraints across λ.
**Catches**: a term nobody decided about; mass or constraint mismatch, which breaks the
exchange validity rather than just the efficiency.

## GATE 1 — λ=1 identity · BLOCKING · seconds

Build the scaled system at λ=1 (any hot selection) and compare per-force-group energies
against the unscaled system on the same frames.

**Pass**: every force group agrees to ≤1e-6 relative.
**Catches**: the scaling code corrupting the topology. The λ=1 replica *is* your answer,
so this is the one gate where failure means a wrong result rather than a slow one.

## GATE 2 — all-hot uniform scaling · BLOCKING · seconds

Mark **every** particle solute — solute, water, ions — and build at λ=0.5. With no cold
region, every term linear in `epsilon` or quadratic in charge must halve and every
invariant term must be unchanged.

**Pass**: each force group's ratio within 1e-6 of 0.5 or 1.0, **and no force group
unaccounted for**.
**Catches**: any term left unscaled — CMAP, pair-specific LJ, exceptions, an exotic
bonded term. This is the gate that generalises across forcefields: it requires no
knowledge of which terms exist, only that they all move together.

Report the full per-group ratio table even on pass. A group that comes back 1.0 where you
expected 0.5 is a decision to make, not a tolerance to loosen.

## GATE 3 — pp/pw/ww partition · BLOCKING · minutes

Gate 2 cannot test the partition, because with everything hot there is no `E_pw`. Measure
the three blocks exactly, using bilinearity, on the **unscaled** system with the
production hot selection. Evaluate only the `NonbondedForce` group:

1. `U_full` — unmodified.
2. `U_a` — solvent charges, `sigma`/`epsilon`, **and all solvent-involving exception
   parameters** zeroed. This leaves `E_pp`.
3. `U_b` — the same for solute. This leaves `E_ww`.
4. `E_pw = U_full − U_a − U_b`, exact, including PME reciprocal space and the dispersion
   correction.
5. Add the scaled bonded groups (solute torsions, CMAP) to `E_pp`.

Then for each production λ compare the measured total potential against
`λ·E_pp + √λ·E_pw + E_ww`, frame by frame.

**Pass**: RMS deviation ≤0.01 kJ/mol per λ.
**Catches**: charges scaled by λ instead of `√λ` — the most common REST2 error, which
puts solute–solvent on λ and solute–solute on λ²; the wrong particle set hot; hot–cold
off-diagonal LJ unscaled.

## GATE 4 — exchange machinery · BLOCKING · hours

**4a, degenerate.** Give every replica the *same* λ, varying only seed and starting
structure. Acceptance must be 1.0.
**Pass**: ≥0.99 in double precision. Any shortfall on GPU must be explained by precision
— neighbouring replicas recompute energies non-reproducibly — not assumed to be.

**4b, non-degenerate.** With the production ladder, compute
`⟨min(1, exp(−Δ))⟩` offline from the actual cross-replica energies and compare to the
observed acceptance per neighbouring pair.
**Pass**: agreement within a few absolute percentage points, pair by pair.
**Catches, and this is the important one**: observed acceptance ≈1.0 across a ladder of
genuinely different Hamiltonians while the prediction is 0.2–0.4 means the criterion is
not seeing your Hamiltonian differences. Configurations are then moving between replicas
without Metropolis weighting and the λ=1 ensemble is biased.

If you use `openmmtools.multistate`, the sampler already computes the full energy matrix,
which makes 4b nearly free and removes this failure class by construction. Verify how it
handles states built on different `System` objects before relying on it.

## GATE 5 — ladder direction · ADVISORY, but a failure invalidates the run's purpose · hours

Compute the order parameter the sampling exists for — R_g, or your slow coordinate — per
replica against effective temperature.

**Pass**: monotonic toward greater disorder.
**Catches**: an anti-correlated ladder. REST2 amplifies solvent cohesion relative to the
solute's effective thermal energy, so compaction at high effective temperature is a known
failure mode for disordered solutes (doi:10.1021/acs.jctc.2c01139). Flat or inverted means
the ladder buys no enhancement for your coordinate while still feeding structures
downward — worse than no ladder.

## GATE 6 — initialization independence · BLOCKING before reporting any ensemble · days

Run two replicates whose replica sets start from deliberately opposite extremes, one from
expanded structures and one from collapsed.

**Pass**: λ=1 estimates agree within combined error, on both the order parameter and the
target observable.
**Catches**: a directionally biased ladder. This is the only gate that detects it. If the
answer depends on where you started, you have not sampled, regardless of how healthy the
exchange statistics look.

## GATE 7 — λ=1 versus plain MD · ADVISORY · days

Independent plain MD at `T_0` with the unscaled system, compared against the λ=1 ensemble
assembled from the demultiplexed replica-index trajectory.

**Pass**: two-sample KS on the potential energy with p > 0.05 at matched effective sample
size, plus agreement on one structural distribution within combined standard error.
Report the autocorrelation time and effective sample size used — a KS test on correlated
frames is meaningless.

---

## Do not conclude convergence from any of these

State explicitly in your report that these were **not** used as evidence:

- Acceptance ratios in the 0.2–0.4 band. They can be ideal while gate 5 or 6 fails.
- An acceptance ratio near 1.0. That is an alarm, not a well-tuned ladder — see gate 4b.
- Round-trip counts. Walkers traverse a ladder that provides no enhancement for your
  coordinate just as happily as one that does.
- Per-replica time-series stability.
- Agreement between seeds that share initial structures. They are fed from the same
  upper rungs.

## Report format

```
GATE 0  System audit              PASS/FAIL   <n forces, n unhandled, n solute particles>
GATE 1  λ=1 identity              PASS/FAIL   <max relative deviation, per-group table>
GATE 2  All-hot scaling           PASS/FAIL   <per-group ratio table>
GATE 3  pp/pw/ww partition        PASS/FAIL   <RMS deviation per λ; Epp/Epw/Eww>
GATE 4a Degenerate acceptance     PASS/FAIL   <acceptance, precision mode>
GATE 4b Predicted vs observed     PASS/FAIL   <per-pair predicted / observed>
GATE 5  Ladder direction          PASS/FAIL   <order parameter vs effective T>
GATE 6  Initialization            PASS/FAIL   <two estimates ± error>
GATE 7  λ=1 vs plain MD           PASS/FAIL   <KS p, n_eff, structural observable>

NOT USED AS EVIDENCE: <list>
UNRESOLVED: <anything you could not determine, and what you tried>
```

Anything you could not verify goes under UNRESOLVED with what you attempted. Do not
infer a gate's result from another gate.
