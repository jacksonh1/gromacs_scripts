# Force fields

How `FF` and `WATER` are resolved, how to add a new force field, and — **read this
before running REST2 with anything new** — how to check that a force field is actually
compatible with the REST2 engine.

Related: [`PARAMETERS.md`](PARAMETERS.md) (all job parameters),
[`../scripts/installation/README.md`](../scripts/installation/README.md) (installing a
force field), [`CATION_PI_WYF.md`](CATION_PI_WYF.md) (cation–π / CHARMM36-WYF — why we
stay on stock `charmm36m`).

---

## Selecting a force field

```bash
FF="amber99sb-ildn"     # default
WATER="tip3p"
```

`FF` is passed to `pdb2gmx -ff`. It is the force-field **directory name minus `.ff`**,
and may be either:

- a force field **bundled with GROMACS** (`amber99sb-ildn`, `charmm27`, …), found in the
  build's own `share/gromacs/top/` (`ls` it to see what your build ships — the 2023.5 and
  2024.3 builds here do **not** include `amber14sb`); or
- one **installed under `GMXLIB`** (exported by `site_config.sh`), found in addition to
  the bundled ones; or
- an **alias** from `FF_ALIASES` in `site_config.sh`.

### Aliases

Ports ship with the release date in the directory name, which is correct but a mouthful.
`site_config.sh` maps a short name to the installed directory:

```bash
declare -A FF_ALIASES=(
  [charmm36m]="charmm36-feb2026_cgenff-5.0"
  [a99sb-disp]="a99SBdisp-25e729d"
)
```

A port with no upstream release name (a99SB-disp is a GitHub repo, not a dated release) is
installed under the **pinned commit** instead — `a99SBdisp-25e729d.ff` — so the run record
still says exactly which parameters were used.

The engine resolves the alias at STEP 1, logs the expansion, and writes the **resolved**
name into `parameters.txt`. The alias is input sugar only — the run record always names
the exact release, so re-pointing an alias at a newer port later cannot silently change
what an old `parameters.txt` means. Both spellings work in a job:

```bash
FF="charmm36m"                       # alias
FF="charmm36-feb2026_cgenff-5.0"     # equivalent
```

**Never rename an installed force-field directory** to a prettier generic name — that
puts the ambiguity back into the run record, which is the thing the alias avoids.

### Water follows the force field

`WATER` is resolved **inside** the force-field directory, so `tip3p` means standard
TIP3P under AMBER and the CHARMM-modified (LJ-on-H) TIP3P under CHARMM, automatically.
Each force field is validated only with its matched water — never cross them.

`WATER` must be a model listed in that force field's `watermodels.dat`; the REST2 engine
checks this at STEP 1 (`scripts/simulation/resolve_water.py`) and fails before any compute
otherwise. The same step handles two things `pdb2gmx`/`solvate` do not do for you:

- **Force-field-specific water names.** `pdb2gmx -water` only accepts GROMACS's built-in
  names (`spc spce tip3p tip4p tip4pew tip5p tips3p`). A model such as a99SB-disp's
  `a99SBdisp_water` is selected through `-water select` with its menu number, and the
  engine then confirms the topology includes `<WATER>.itp`.
- **The solvent box must match the model's site count.** A 4-site water in the 3-site
  `spc216.gro` box is a fatal grompp mismatch. The box is the force field's own
  `<WATER>.gro` if it ships one (its atom names match the `.itp`), else `spc216.gro` /
  `tip4p.gro` / `tip5p.gro` by site count. 3-site built-ins resolve to exactly the old
  `-water <WATER>` + `spc216.gro`, so existing TIP3P runs build byte-identically.

### CHARMM changes the mdp

Any `FF` starting with `charmm` switches the generated mdp to force-switched van der
Waals and forces the cutoff:

```
vdwtype      = cutoff
vdw-modifier = force-switch
rvdw-switch  = 1.0
rvdw = rcoulomb = 1.2      (CUTOFF_NM forced to 1.2, with an [INFO] notice)
DispCorr     = no
```

CHARMM36/36m was parameterized this way; running it with the AMBER plain-cutoff settings
runs fine and looks plausible but gives subtly wrong forces. The AMBER branch is
unchanged and byte-identical to before this gating existed.

**If you add a force field that is neither AMBER-like nor CHARMM-like, check its
published nonbonded settings first** — the mdp generator only knows these two cases, and
will silently give a new force field the AMBER treatment.

---

## Force fields and REST2 — check before you run

**T-REMD and plain MD accept any force field `pdb2gmx` can build.** They use the
topology exactly as written and scale nothing, so there is nothing to verify.

**REST2 is different, and the failure is silent.** REST2 samples by scaling the
*solute's* Hamiltonian by λ, which the pipeline does by having
`plumed partial_tempering` rewrite the processed topology. That script understands a
fixed set of topology sections and **passes everything else through untouched, without
warning**. A force field whose solute energy includes a term outside that set is
therefore only *partially* scaled: the run completes, exchanges at plausible rates, and
samples an ensemble that is not the one it reports.

### What `partial_tempering` scales

Measured by diffing a λ=1.0 topology against a λ=0.5 one:

| Section | Scaled? | What it is |
|---|---|---|
| `[ atoms ]` | **yes** | per-atom charges (×√λ) |
| `[ atomtypes ]` | **yes** | LJ ε of the marked types |
| `[ nonbond_params ]` | **yes** | explicit LJ pair parameters |
| `[ pairtypes ]` | **yes** | 1-4 pair parameters |
| `[ dihedrals ]` | **yes** | torsion barriers |
| `[ bonds ]` `[ angles ]` | no — **correct** | stiff terms, deliberately unscaled in REST2 |
| `[ dihedraltypes ]` | no — **correct** | commented out; the scaled values are inlined into `[ dihedrals ]` |
| `[ settles ]` `[ exclusions ]` | no — **correct** | water constraints / bookkeeping, no solute energy |
| **anything else** | **no — and that is the hazard** | e.g. `[ cmap ]` |

Bonds and angles being unscaled is the method, not a bug: REST2 scales charges, LJ and
torsions, because those govern conformational sampling. The danger is only a term that
*does* govern conformation and is not in the scaled list.

### The rule

> A force field is REST2-compatible **iff every solute-energy term it uses that governs
> conformation is one of: charges, LJ, 1-4 pairs, or dihedrals.** Any additional
> conformational term — CMAP, tabulated torsions, polarization — is passed through
> unscaled and breaks the method.

### The known failure: CHARMM CMAP

CHARMM's CMAP is a backbone φ/ψ cross-term correction. It appears as `[ cmap ]` and
`[ cmaptypes ]`, `partial_tempering` has no handling for it at all, and it applies
**only to protein** — i.e. exactly the hot region. Measured on a CHARMM36m topology:
`[ cmap ]` (41 entries) and `[ cmaptypes ]` (1475 entries) are **byte-identical between
λ=1.0 and λ=0.5**.

So a CHARMM REST2 run would scale charges, LJ and dihedrals while leaving the backbone
conformational term at full strength. **Neither existing self-check catches it:** λ=1.0
is still exact, so the `scale=1.0 → P=1.0` sanity pair passes; exchanges still happen at
normal rates, so the acceptance gate passes.

The REST2 engine therefore **rejects `FF=charmm*` at STEP 1**. This is a limitation of
the scaling tool, not of the GROMACS build (the 2023.5 REST2 build reads `GMXLIB` and
builds CHARMM topologies fine) and not of REST2 as a method.

### The automatic check (STEP 7)

The name check above only catches CHARMM, and only because someone knew to write it. So
every REST2 run also checks the scaled topologies **by content** at STEP 7:
`scripts/simulation/check_rest2_topology.py` compares `topol_rep000.top` (λ=1) with the
lowest-λ replica and fails the job unless

- every hot charge is ×√λ, every hot atom type's LJ ε is ×λ, and every
  `[ nonbond_params ]` / `[ pairtypes ]` override is ×λ (both types hot) or ×√λ (one hot);
- every torsion's force constants are ×√λ per hot atom in positions 1 and 4;
- bonds, angles and everything else are unchanged; and
- no section `partial_tempering` cannot scale (`[ cmap ]`, `[ cmaptypes ]`,
  polarization, restraints, …) and no `[ dihedrals ]` function type outside
  {1,2,3,4,5,9} is present.

It checks numbers, not just "did the section change", so a partly applied scaling fails as
well as a missing one. On CHARMM36m it fails on exactly `[ cmaptypes ]` and `[ cmap ]`;
on amber99sb-ildn and a99SB-disp it passes. Its output is kept in
`topol/check_rest2_topology.log`.

### Known status

| Force field | REST2 | Why |
|---|---|---|
| `amber99sb-ildn` | **supported** | no CMAP; all solute terms are scaled |
| `a99sb-disp` (`a99SBdisp-25e729d`), `WATER=a99SBdisp_water` | **supported — validated 2026-09-28** | no CMAP; backbone O–H `[ nonbond_params ]` override scaled correctly; see the validation record below |
| `amber14sb` | supported in principle (no CMAP) | not installed in the builds here; run the check below if you install a port |
| `charmm36m` (any `charmm*`) | **rejected at STEP 1** (and STEP 7) | CMAP, unscaled |
| AMBER **ff19SB** | **rejected at STEP 7** | CMAP, unscaled; the cluster's `amber19sb.ff` also needs GROMACS ≥2026 syntax the REST2 build cannot read |
| Drude / polarizable | **rejected at STEP 7** | polarization terms are not scaled |
| anything else new | **unknown — validate it (below)** | |

### Validating a new force field for REST2

Two independent questions, both once per force field, both CPU-only:

1. **Is the scaling right?** `scripts/validation/validate_rest2_scaling.sh` builds the
   λ=1, all-atoms-hot and protein-hot topologies from a solvated system and checks every
   energy term's ratio against the unscaled topology, frame by frame, using energy groups
   for the protein/water split:

   ```bash
   bash scripts/validation/validate_rest2_scaling.sh em.gro system.top /tmp/v2 0.5 [traj.xtc]
   ```

   | check | expected |
   |---|---|
   | (a) protein hot, λ=1 | every term ×1 (tol 1e-6) |
   | (b) all atoms hot, λ | bonds/angles ×1, every other term (incl. dispersion correction, PME reciprocal) ×λ |
   | (c) protein hot, λ | protein–protein ×λ, protein–water ×√λ, water–water ×1, torsions and 1-4 ×λ |

2. **Is the port itself right?** A GROMACS port is a hand conversion, and a single wrong
   torsion goes unnoticed by (1). Compare it term by term against an independent
   implementation of the same force field on identical coordinates (for AMBER-family
   force fields: tleap + sander from AmberTools), over tens of varied conformations, and
   compare the per-torsion parameters to find the exact line when a term disagrees. See
   the a99SB-disp record below for what that looks like and what to expect.

The manual section-diff below is what the STEP 7 tool automates; keep it for looking at
an unfamiliar topology by eye.

### Manual section-diff (reference)

It takes a couple of minutes and needs no GPU.

```bash
REPO=/path/to/gromacs_REMD
NEWFF="your-force-field-name"          # as passed to pdb2gmx -ff
PDB="$REPO/example/input_pdbs/helix_fusion.pdb"

mkdir -p /tmp/ffcheck && cd /tmp/ffcheck
source "$REPO/site_config.sh"
export LD_LIBRARY_PATH="${HCOLL_COMPAT_DIR}:${LD_LIBRARY_PATH:-}"
source "$REST2_GMXRC"                   # the REST2 build, not the T-REMD one
source "$PLUMED_SH"

# 1. build a topology and a processed (fully expanded) one
gmx_mpi pdb2gmx -f "$PDB" -o c.gro -p c.top -i p.itp -ff "$NEWFF" -water tip3p -ignh
gmx_mpi editconf -f c.gro -o box.gro -bt dodecahedron -d 1.0
printf 'integrator = md\ndt = 0.002\nnsteps = 0\ncutoff-scheme = Verlet\nnstlist = 10\ncoulombtype = PME\nrcoulomb = 1.2\nrvdw = 1.2\npbc = xyz\n' > pp.mdp
gmx_mpi grompp -f pp.mdp -c box.gro -p c.top -pp processed.top -o pp.tpr -maxwarn 5

# 2. mark the solute and scale it at two lambdas
python3 "$REPO/scripts/simulation/mark_hot_region.py" processed.top > marked.top
plumed partial_tempering 1.0 < marked.top > pt_1.0.top
plumed partial_tempering 0.5 < marked.top > pt_0.5.top

# 3. which sections actually changed?
python3 - <<'EOF'
def secmap(path):
    cur=None; out={}
    for i,l in enumerate(open(path)):
        s=l.strip()
        if s.startswith('[') and s.endswith(']'): cur=s
        out[i]=cur
    return out
a=open('pt_1.0.top').read().splitlines()
b=open('pt_0.5.top').read().splitlines()
assert len(a)==len(b), "topologies differ in length — investigate before trusting this"
secs=secmap('pt_1.0.top')
from collections import Counter
changed=Counter(secs[i] for i,(x,y) in enumerate(zip(a,b)) if x!=y)
present={v for v in secs.values() if v}
print("CHANGED (scaled):")
for k,v in changed.most_common(): print(f"  {v:8d}  {k}")
print("\nUNCHANGED sections present in the topology:")
for k in sorted(present - set(changed)): print(f"            {k}")
EOF
```

**Reading the result.** The `CHANGED` list must include `[ atoms ]`, `[ atomtypes ]`,
`[ dihedrals ]`, and whichever of `[ nonbond_params ]` / `[ pairtypes ]` the force field
uses. Then go through the `UNCHANGED` list and ask of each: *does this section carry
solute energy that affects conformation?*

- `[ bonds ]` `[ angles ]` `[ bondtypes ]` `[ angletypes ]` — expected, fine.
- `[ dihedraltypes ]` — expected, fine. `partial_tempering` comments this block out and
  writes the scaled torsion parameters **inline** into `[ dihedrals ]`, so the type
  block is identical across λ while the actual barriers do scale. Confirm by checking
  that `[ dihedrals ]` is in the CHANGED list; if `[ dihedraltypes ]` is unchanged *and*
  `[ dihedrals ]` did not change, nothing was scaled and something is wrong.
- `[ settles ]` `[ exclusions ]` `[ defaults ]` `[ moleculetype ]` `[ system ]`
  `[ molecules ]` `[ pairs ]` — bookkeeping or water, fine.
- **Anything else is a red flag.** `[ cmap ]`, `[ cmaptypes ]`, `[ polarization ]`,
  `[ thole_polarization ]`, `[ pairs_nb ]`, tabulated dihedral types — these carry
  energy and are not being scaled. Do not run REST2 with that force field.

For reference, a CHARMM36m topology produces exactly this red flag:

```
CHANGED (scaled):
    151434  [ pairtypes ]
       865  [ dihedrals ]
       565  [ atomtypes ]
       433  [ atoms ]
       333  [ nonbond_params ]

UNCHANGED sections present in the topology:
            [ angles ]
            [ angletypes ]
            [ bonds ]
            [ bondtypes ]
            [ cmap ]            <-- RED FLAG
            [ cmaptypes ]       <-- RED FLAG
            [ defaults ]
            [ dihedraltypes ]   (fine — scaled torsions are inlined into [ dihedrals ])
            [ exclusions ]
            [ molecules ]
            [ moleculetype ]
            [ pairs ]
            [ settles ]
            [ system ]
```

### If the check fails

You have three options, in order of effort:

1. **Use T-REMD instead.** It scales nothing, so every force field is handled exactly as
   `pdb2gmx` wrote it. This is the right answer almost every time.
2. **Pick a force field that passes.** AMBER without CMAP.
3. **Teach `partial_tempering` to scale the missing term**, then re-validate — including
   a fresh λ=1.0 sanity check and an acceptance-rate check. Only worth it if the science
   genuinely requires both that force field and solute tempering.

### If you add a supported force field

Add it to the "Known status" table above and to `FF_ALIASES` in `site_config.sh`, and
record what you measured. The next person should not have to rediscover it.

---

## Validation record: a99SB-disp (2026-09-28)

**Port:** `github.com/paulrobustelli/Force-Fields`, `Gromacs_FFs/a99SBdisp.ff` at commit
`25e729d` (2025-01-27), installed unmodified as `$GMXLIB/a99SBdisp-25e729d.ff`. The
authors' definitive version is the Desmond/viparr one in the same repo; they state the
GROMACS port is "not guaranteed bug-free", and its history has real fixes (NTHR angles
2021, CPRO improper 2022, HIP improper 2024, ion LJ Oct 2024). Use it with its own water,
`WATER=a99SBdisp_water` (TIP4P-D with increased dispersion; 4-site) and the AMBER
nonbonded settings (1.0 nm, PME, `DispCorr=EnerPres` — the published non-Desmond runs
used a 10 Å cutoff with PME under Amber's default long-range dispersion correction).

**Port fidelity** — GROMACS 2023.5 `mdrun -rerun` vs AmberTools `sander` (Amber-format
a99SB-disp from the same repo), identical coordinates, vacuum, no cutoff, 51 frames from
500 K MD each:

| term | 23-residue peptide (all 20 aa + HID/HIE/HIP, termini) | gp130 binder (2 disulfides) |
|---|---|---|
| bond | 2.4e-5 | 3.7e-5 |
| angle | 8.3e-6 | 7.1e-6 |
| torsions (proper + improper) | 7.7e-6 | 8.5e-6 * |
| LJ-14 | 8.3e-6 | 6.3e-6 |
| Coulomb-14 | 1.4e-5 | 3.8e-5 |
| LJ | 1.2e-5 | 1.2e-5 |
| Coulomb | 5.7e-5 | 4.8e-5 |

(max relative difference over all frames; the ~1e-5 floor is single precision plus the
two codes' Coulomb constants differing by 3.5e-5 — the amber99sb-ildn control gives the
same floor.) Findings worth knowing:

- \* **tleap, not the port, was wrong on serine.** On gp130, sander's torsion energy was
  up to 4 kJ/mol off. Per-torsion comparison put it on the serine `OG-CB-CA-C` torsion:
  tleap had added the generic `X-CT-C8-X` 3-fold term on 8 of 9 serines as well as the
  specific `C-CT-C8-OH` 1-fold term. The Desmond definition has only the specific term
  (viparr: exact match beats wildcard), and so does the GROMACS port. With tleap's extra
  term removed the gp130 torsions agree to 8.5e-6. **If you validate an AMBER-family port
  against tleap, judge disagreements against the Desmond/original definition.**
- **The backbone O–H override applies to 1-4 pairs in both codes.** Amber's `LJEDIT
  OB HB` and GROMACS's `[ nonbond_params ] OB HB` both enter the generated 1-4 LJ
  (removing either shifts LJ-14 by the same 0.036 kJ/mol on the peptide), so O(i)–H(i+1)
  is treated identically. Under REST2 the scaled `OB_ HB_` entry therefore scales those
  1-4 terms too, as it should.
- **HIP works** on this build despite upstream issue #25; the HIP improper is listed in a
  different atom order than tleap's but is energetically equivalent (within 0.01 kJ/mol
  over the hot frames).
- **Water** matches the Desmond definition (charges H +0.59 / M −1.18, O σ 3.165 Å
  ε 0.238764 kcal/mol, M-site 0.131937768, 0.9572 Å / 104.52°).
- **Control:** GROMACS's *stock* `amber99sb-ildn` vs tleap `ff99SBildn` shows the same
  ~1e-5 floor on every term except torsions, which carry a **constant** −74.64 kJ/mol
  offset (identical in every frame). That is a representation difference in the stock
  port's torsion constants — no force, no exchange-criterion effect — not an error.

**REST2 scaling** (`validate_rest2_scaling.sh`, solvated gp130 in `a99SBdisp_water`,
λ=0.5): (a) λ=1 equals the unscaled topology on all 22 terms to 1e-6; (b) all-hot: bonds
and angles 1.000000, everything else 0.500000–0.500001; (c) protein-hot: protein–protein
0.500000, protein–water 0.707107, water–water 1.000000, torsions and 1-4 0.500000. The
same script passes on amber99sb-ildn/TIP3P at λ=0.35.

**Engine runs** (gp130, 300–600 K effective ladder, 24 replicas, 1 L40S, NPT, 2 ns; data in
`data/rest2_a99sbdisp_validation/`):

| | amber99sb-ildn / TIP3P (pilot `rep24_gpu1`) | a99SB-disp / `a99SBdisp_water` |
|---|---|---|
| mean neighbour acceptance (min–max, 23 pairs) | 0.449 (0.38–0.50) | 0.443 (0.39–0.51) |
| ns/day per replica | 162 | 114 (4-site water, ~30% slower) |
| λ=1 backbone RMSD mean / max | 2.6 / 4.1 Å | 1.7 / 2.4 Å |

STEP 7's `check_rest2_topology.py` passed; no LINCS/SETTLE warnings; λ=1 replica at
300.0 K. System density 1025 kg/m³ (vs 1013 with TIP3P, which is known to be under-dense);
**pure `a99SBdisp_water`: 995.1 ± 0.7 kg/m³ at 300 K, 1 bar** (experiment 996.5).

**Identical-λ exchange check** (2 replicas, both λ=1): acceptance 84% on GPU, not 100%.
This is GPU energy noise in PLUMED's hrex re-evaluation (`dplumed` sd 0.36 kT), identical
for amber99sb-ildn (80%) and absent on CPU (exactly 1.0) — see `knowledgebase/GOTCHAS.md`,
"on GPU, hrex exchange energies carry ~0.37 kT of random noise". It is not specific to
a99SB-disp and applies to every GPU REST2 run.
