# Force fields

How `FF` and `WATER` are resolved, how to add a new force field, and — **read this
before running REST2 with anything new** — how to check that a force field is actually
compatible with the REST2 engine.

Related: [`PARAMETERS.md`](PARAMETERS.md) (all job parameters),
[`../scripts/installation/README.md`](../scripts/installation/README.md) (installing a
force field).

---

## Selecting a force field

```bash
FF="amber99sb-ildn"     # default
WATER="tip3p"
```

`FF` is passed to `pdb2gmx -ff`. It is the force-field **directory name minus `.ff`**,
and may be either:

- a force field **bundled with GROMACS** (`amber99sb-ildn`, `amber14sb`, `charmm27`, …),
  found in the build's own `share/gromacs/top/`; or
- one **installed under `GMXLIB`** (exported by `site_config.sh`), found in addition to
  the bundled ones; or
- an **alias** from `FF_ALIASES` in `site_config.sh`.

### Aliases

Ports ship with the release date in the directory name, which is correct but a mouthful.
`site_config.sh` maps a short name to the installed directory:

```bash
declare -A FF_ALIASES=(
  [charmm36m]="charmm36-feb2026_cgenff-5.0"
)
```

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

### Known status

| Force field | REST2 | Why |
|---|---|---|
| `amber99sb-ildn`, `amber14sb` | **supported** | no CMAP; all solute terms are scaled |
| `charmm36m` (any `charmm*`) | **rejected at STEP 1** | CMAP, unscaled |
| AMBER **ff19SB** | **do not use** | also uses CMAP — same failure, **not** caught by the `charmm*` guard |
| Drude / polarizable | **do not use** | polarization terms are not scaled |
| anything else new | **unknown — run the check below** | |

### How to check a new force field

Do this **once per force field**, before any production REST2 run. It takes a couple of
minutes and needs no GPU.

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
