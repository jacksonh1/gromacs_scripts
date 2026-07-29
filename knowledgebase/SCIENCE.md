# Science

Rationale for the simulation protocols. Terse. Operational flags/pitfalls are in
`CLAUDE.md`; this is the *what and why* of the physics.

---

## What these pipelines are for

Every run starts from a **folded/designed input pose — never an unfolded or
extended state.** The input pose is the *reference* the analysis is measured
against (RMSD = drift from the design, RMSF = local flexibility). Four objectives:

1. Stability/rigidity of a given starting structure
2. Identifying flexible regions (per-residue)
3. Variant comparison — same protocol across variants; which best retains its pose
4. Bound-state sampling — a complex in its bound pose

Engine mapping: plain **MD** → mainly #4 (+ optionally #1, #2); **T-REMD** and
**REST2** → primarily #1–3. Neither is a folding-from-unfolded tool.

## T-REMD (temperature replica exchange)

N replicas span a geometric temperature ladder (`T_MIN`…`T_MAX`). Coordinates
exchange between adjacent temperature slots via a Metropolis criterion; the high
temperatures cross barriers, the exchanges feed that enhanced sampling back down
to the temperature of interest.

**Key concept — slots, not configurations.** Each `prod/rep{i}/` is a *fixed
temperature slot*. `rep000` is the constant-`T_MIN` (e.g. 300 K) thermodynamic
ensemble — analyze it directly. Under NPT, GROMACS adds the *PV* term to the
exchange criterion automatically.

## REST2 (replica exchange with solute tempering) — in progress

Instead of heating the whole box, REST2 runs **every replica at the same physical
temperature** and scales the Hamiltonian of a chosen *solute* region per replica.
Here the solute is the **whole protein**; water/ions are unscaled.

Scaling by λ_i = `T_MIN`/T_i (geometric effective-temperature ladder), applied to
the solute:
- charges × √λ, LJ ε × λ, proper dihedrals × λ (⇒ solute–solute × λ,
  solute–solvent × √λ, solvent–solvent unchanged)

**rep000 has λ=1 → the true, unscaled physical ensemble** (same role as 300 K
rep000 in T-REMD). Because only solute DOF are "heated", REST2 needs **far fewer
replicas** than T-REMD for the same effective range (~8–16 vs ~48). Implemented
with the PLUMED-patched GROMACS: `plumed partial_tempering` builds the per-replica
scaled topologies, `mdrun -hrex` runs Hamiltonian exchange. See
`plans/REST2-pipeline.md`.

## Plain MD

Single-temperature NPT production (V-rescale thermostat + C-rescale barostat).
Position restraints (`-DPOSRES`) held through EM/heat/density to preserve the
input pose, released at production start (or after an optional unrestrained
`relax/` stage if `RELAX_NS>0`). Mainly bound-state equilibrium sampling.

## Ensembles

- **REMD production:** NVT (default, constant volume) or NPT (`ENSEMBLE=NPT`,
  C-rescale barostat). The constant-temperature-slot interpretation is unaffected.
- **Plain MD production:** always NPT.
