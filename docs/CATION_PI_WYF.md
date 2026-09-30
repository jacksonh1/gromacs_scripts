# Cation–π accuracy and the CHARMM36-WYF question

Whether we can use the CHARMM36 cation–π correction (commonly called **WYF**) with our
installed `charmm36m` port, what it would take, and why — as of 2026-09-30 — the
recommendation is **not to attempt it**.

Related: [`FORCE_FIELDS.md`](FORCE_FIELDS.md) (how `FF`/`WATER` resolve, adding a force
field, the CHARMM mdp gate).

---

## TL;DR

- Having `charmm36m` installed does **not** give us WYF. WYF is an **add-on NBFIX patch**,
  not part of any stock CHARMM36 release.
- Our port (`charmm36-feb2026_cgenff-5.0.ff`) ships NBFIX for halogens, metals, lipids,
  carbohydrates, ethers and CGenFF — but **no cation–π (cation–aromatic) terms**. Verified
  by scanning its `nbfix.itp`.
- **No authoritative ready-made GROMACS-format WYF port was found.** MacKerell distributes
  WYF only in **CHARMM format** (toppar). Adding it to GROMACS means hand-transcribing
  NBFIX lines and mapping atom types — a silent-wrong-results risk.
- **Recommendation: stay on stock `charmm36m`.** Only revisit if an authoritative,
  atom-type-matched GROMACS WYF port appears *and* we validate it against a published
  cation–π number first.

---

## Background: what WYF is

Standard fixed-charge force fields (CHARMM36 included) systematically **underestimate
cation–π interactions** — the attraction between a cation (e.g. a protonated Lys/Arg
sidechain, a methylated ammonium, a Na⁺) and the face of an aromatic ring (Trp/Tyr/Phe).

The fix is a set of **NBFIX** terms. NBFIX ("nonbonded fix") overrides the Lennard-Jones
`sigma`/`epsilon` for **specific atom-type pairs**, replacing the value the standard
combination rule would produce. WYF applies NBFIX between cation atom types and the
**aromatic ring carbons** of **W**/**Y**/**F** — hence the name. It touches nonbonded
parameters only; no bonded terms, no CMAP, and (in principle) no mdp change.

Key reference: Khan, MacKerell & Reuter, *Cation-π Interactions between Methylated
Ammonium Groups and Tryptophan in the CHARMM36 Additive Force Field*, JCTC 2019.

---

## MD primer: the "12-6 cutoff" and "force-switch" (for non-MD-experts)

This matters because a cation–π NBFIX is only correct under the **same van der Waals
cutoff scheme it was fitted with**.

- **The 12-6 / Lennard-Jones potential.** Two non-bonded atoms feel a van der Waals
  interaction: a steep repulsion `~1/r¹²` (atoms can't overlap) plus a weak attraction
  `~1/r⁶` (stickiness). "12-6" = those exponents. Each atom pair is described by two
  numbers: `sigma` (size) and `epsilon` (well depth). NBFIX overrides exactly these two
  for chosen pairs.

- **The cutoff.** vdW is tiny but never exactly zero at distance. Computing every pair to
  infinity is too expensive, so it is cut off at a radius (here 1.2 nm) and ignored beyond.
  Naively chopping the potential there leaves a **kink** — a discontinuous force at the
  cutoff — which causes energy drift and artifacts.

- **The cutoff "scheme" = how that kink is smoothed:**
  - *Plain cutoff:* chop (or simply shift) at the radius. Common in AMBER-style setups.
  - *Force-switch:* over a window (e.g. 1.0 → 1.2 nm) smoothly ramp the force to exactly
    zero at the cutoff, so the force stays continuous. **CHARMM force fields were fitted
    with force-switch**, and their parameters assume it. Running CHARMM with a plain cutoff
    is silently wrong — which is why this repo's mdp **auto-gates on `FF==charmm*`** and
    forces `CUTOFF_NM=1.2` with force-switch vdW (see `FORCE_FIELDS.md` and the CHARMM
    gotcha in `../knowledgebase/GOTCHAS.md`).

- **Why it matters for WYF.** The NBFIX `sigma`/`epsilon` overrides were fitted under some
  specific cutoff scheme. Applying them under a *different* scheme changes the effective
  interaction from what the authors calibrated — **even with perfectly correct atom
  types**. So the scheme WYF assumes must be confirmed before trusting any numbers. As of
  this writing that confirmation was **not** obtained (the paper and MacKerell page were
  not retrievable without access).

---

## Research findings (2026-09-30)

1. **Our port has no cation–π NBFIX.** `charmm36-feb2026_cgenff-5.0.ff/nbfix.itp` contains
   only halogen (CLGR1/BRGR1), CHARMM-metal, lipid, carbohydrate, ether and CGenFF terms.
   No cation↔aromatic-carbon pairs.

2. **No stock GROMACS port bundles WYF.** The same absence was confirmed in the widely used
   `intbio/charmm36-mar2019` port, which is also what the MacKerell lab distributes as its
   GROMACS port. WYF is distributed by MacKerell **in CHARMM format only**.

3. **The community treats this as unsolved/fiddly.** An official GROMACS forum thread
   (Aug 2023) asks exactly how to add WYF to a GROMACS port and raises the same concerns we
   have — which file to edit (`nbfix.itp` vs `ffnonbonded.itp`) and whether the atom-type
   pairs even exist — with no clean, authoritative answer.

4. **Provenance is tangled.** "WYF" refers to the Trp/Tyr/Phe cation–π NBFIX set; the
   Khan/MacKerell/Reuter 2019 paper covers the methylated-ammonium–Trp piece specifically.
   The full WYF set may combine sources, which makes "download the one correct file" harder,
   not easier.

---

## Risk assessment (against a "zero-risk" bar)

Two concerns were raised, ranked above having the feature at all:

- **Concern A — don't break the existing `charmm36m` FF.** *Fully avoidable.* Any attempt
  would `cp -r` the port to a new directory, edit only the copy, add a **separate** alias
  (e.g. `charmm36m-wyf`), and `sha256sum` every file of the original before/after to prove
  it is byte-identical. The original alias keeps resolving to an unchanged directory;
  existing and in-flight jobs see nothing different. **Not the blocker.**

- **Concern B — download the right files and patch correctly.** *Cannot be made
  risk-free here.* There is no authoritative GROMACS artifact to drop in; it requires
  manual CHARMM→GROMACS transcription plus atom-type mapping. Failure modes:
  - Wrong/undefined atom type → `grompp` **fatal error** (loud, safe).
  - Atom type exists but maps to the **wrong pair** → applied **silently wrong** (dangerous).
  Add the unverified vdW-scheme assumption and the tangled provenance, and the
  silent-wrong surface exceeds a zero-risk tolerance. It also conflicts with this repo's
  "fail loudly / no silent wrong results" philosophy.

---

## Would a different force field be easier? (No)

A natural follow-up: rather than patch WYF, is there a force field — an AMBER one, say —
that handles cation–π better *and* is easier to use in GROMACS? Checked 2026-09-30; the
answer is **no easy win**.

- **Cation–π is underestimated by *every* fixed-charge force field** — AMBER (ff14SB,
  ff19SB), CHARMM, OPLS alike. It is not a CHARMM-specific gap. The physical reason is that
  most of the cation–π attraction comes from **polarization/induction**, which fixed-charge
  models structurally lack.

- **AMBER does not help on this axis.** ff14SB/ff19SB now have clean, validated GROMACS
  ports (AMBER14SB/AMBER19SB), and AMBER is arguably *easier to run* in GROMACS than CHARMM
  (plain cutoff, no force-switch gate, simpler mdp). But its cation–π is **just as
  underestimated**, and there is **no clean cation–π-specific AMBER drop-in** equivalent to
  WYF. The well-known AMBER NBFIX work (**CUFIX**, Aksimentiev lab) targets **charged-group
  over-attraction** (amine–carboxylate, amine–phosphate, ion pairing) — *not* cation–π.

- **The genuine fixes are polarizable force fields** — CHARMM **Drude-2013** (improved
  cation–π, Lin/Roux 2020) and **AMOEBA**. These are physically correct for cation–π, but
  they are **harder, not easier** in GROMACS: Drude support is limited/fragile and AMOEBA is
  effectively Tinker/OpenMM, not GROMACS. They fail the "easier to use" bar decisively.

| Option | cation–π | GROMACS ease |
| --- | --- | --- |
| stock `charmm36m` / AMBER (ff14SB, ff19SB) | underestimates | easy |
| CHARMM36-WYF NBFIX | better | hard (hand-patch — see above, rejected) |
| CUFIX (AMBER/CHARMM) | not cation–π (wrong target) | medium |
| Drude / AMOEBA (polarizable) | correct | very hard / unsupported in GROMACS |

**Practical caveat for our use.** Because any fixed-charge force field under-weights
cation–π, when cation–π contacts matter to a design's binding, treat them in our MD
**qualitatively, not quantitatively** — in particular for bound-state sampling (objective
#4). A weak or transient cation–π contact in a trajectory may be real but under-represented;
do not read its occupancy as a calibrated interaction strength.

---

## Recommendation

**Stay on stock `charmm36m`.** Do not hand-patch WYF in. No alternative force field offers
better cation–π *and* easier GROMACS use (see the table above).

Revisit only if **both** of the following are satisfied:

1. An **authoritative, atom-type-matched GROMACS-format WYF port** is obtained — either a
   vetted published port, or files supplied directly by the MacKerell lab — with the
   intended **vdW cutoff scheme stated**; and
2. It is **validated** by reproducing a cation–π quantity reported in the source (e.g. an
   interaction free energy / PMF) before any production run.

If we do reach that point, install it as a *new* aliased port (`charmm36m-wyf`) per
Concern A above — never by editing the existing directory.

### Open item (tabled)

Confirming the exact modified atom-type list and the vdW cutoff scheme WYF assumes
requires the paper / MacKerell toppar, which were not retrievable in this session. The user
has paper access; this was **deliberately tabled**. Pick it up here if the question is
reopened.

---

## Sources

- GROMACS forum — *CHARMM36m WYF parameters*:
  https://gromacs.bioexcel.eu/t/charmm36m-wyf-parameters/6957
- Khan, MacKerell & Reuter 2019 (cation–π, JCTC):
  https://www.researchgate.net/publication/329768707_Cation-p_Interactions_between_Methylated_Ammonium_Groups_and_Tryptophan_in_the_CHARMM36_Additive_Force_Field
- `intbio/gromacs_ff` charmm36-mar2019 port (no cation–π NBFIX):
  https://github.com/intbio/gromacs_ff/tree/master/charmm36-mar2019.ff
- MacKerell lab force-field page:
  https://mackerell.umaryland.edu/charmm_ff.shtml
- GROMACS CHARMM36 recommended settings (force-switch, rvdw 1.2, rvdw-switch 1.0):
  https://manual.gromacs.org/current/user-guide/force-fields.html
- CUFIX — NBFIX for CHARMM & AMBER, charged-group corrections (Aksimentiev lab):
  https://bionano.physics.illinois.edu/CUFIX
- Porting AMBER ff14SB/ff19SB to GROMACS (validated ports):
  https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15006112/v1
- Drude polarizable cation–π improvement (Lin et al. 2020, J. Comput. Chem.):
  https://onlinelibrary.wiley.com/doi/10.1002/jcc.26067
- Polarizable force fields review (arXiv):
  https://arxiv.org/pdf/1910.14237
