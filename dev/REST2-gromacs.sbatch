#!/bin/bash
#SBATCH --job-name=rest2
#SBATCH -p pi_keating
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --gres=gpu:l40s:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

###############################################################################
# GROMACS REST2 (Replica Exchange with Solute Tempering 2) — single-node
# ---------------------------------------------------------------------------
# REST2 scales only the protein (solute) Hamiltonian across replicas while
# keeping all replicas at the same physical temperature T_REF.  Compared to
# T-REMD, this is far more efficient because the many water degrees of freedom
# are not heated — the same effective temperature range requires ~4–5x fewer
# replicas.
#
# Effective Hamiltonian for replica m:
#   H_m(x) = λ_m · E_pp(x)  +  √λ_m · E_pw(x)  +  E_ww(x)
#   λ_m = T_REF / T_eff_m
#
# GROMACS implementation (native free-energy / -hrex framework):
#   A-state : original force-field parameters   (λ=0, reference replica)
#   B-state : ε_B=0, q_B=0 for all protein atoms
#   vdw-lambdas_i  = 1 − T_REF/T_eff_i      → protein-protein ε scales as λ·ε_A
#   coul-lambdas_i = 1 − √(T_REF/T_eff_i)   → protein charges scale so q·q ∝ λ·q²_A
#
# The √λ protein–water LJ cross-term arises automatically from the geometric
# Lorentz–Berthelot mixing rule:
#   ε_cross(λ) = √(ε_protein(λ) · ε_water) = √(1−λ_vdw) · ε_cross_A = √λ_REST2 · ε_cross_A ✓
#
# Note: only non-bonded interactions are scaled (standard practice for
# peptide REST2).  Bonded terms are kept at full strength for all replicas.
#
# Usage:
#   sbatch REST2-gromacs.sbatch
#   sbatch REST2-gromacs.sbatch config.sh      # source a config file first
#   sbatch -n 24 --export=ALL,T_EFF_MAX=400,TOTAL_NS=50,PDB_IN=my.pdb REST2-gromacs.sbatch
###############################################################################

set -euo pipefail

CURRENT_TIME=$(date +%s)
CURRENT_STEP=0
PRESERVE_FROM_STEP=${PRESERVE_FROM_STEP:-9}

# ── Source optional config file ──────────────────────────────────────────────
if [[ -n "${1:-}" && -f "$1" ]]; then
  echo "[INFO] Sourcing config: $1"
  # shellcheck disable=SC1090
  source "$1"
fi

###############################################################################
# STEP 0 — Environment
###############################################################################
echo "===== JOB START $(date) ====="

module purge
module load gcc
module load cuda/12.9.1
module load openmpi/4.1.4

set +u
source "/home/jhalpin/opt/gromacs/2024.2/bin/GMXRC"
set -u
GMX="${GMX:-gmx_mpi}"

command -v "$GMX" >/dev/null 2>&1 || { echo "[ERROR] $GMX not found"; exit 1; }
echo "[OK] GROMACS found: $($GMX --version 2>&1 | head -1)"

###############################################################################
# STEP 1 — Parameters
###############################################################################
CURRENT_STEP=1
echo "===== STEP 1: Setting parameters ====="

PDB_IN="${PDB_IN:-IL7-2-V3-cys.pdb}"
OUTBASE="${OUTBASE:-IL7-2-V3-cys}"
OUTDIR="${OUTDIR:-rest2_${OUTBASE}}"

# Force field & water model
FF="${FF:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
BOX_SHAPE="${BOX_SHAPE:-dodecahedron}"
BOX_BUFFER="${BOX_BUFFER:-1.0}"
NEUTRALIZE="${NEUTRALIZE:-1}"
SALT_MOLAR="${SALT_MOLAR:-0.15}"

# Simulation parameters
REPLICAS=${SLURM_NTASKS:-${REPLICAS:-48}}
DT_PS="${DT_PS:-0.002}"
T_REF="${T_REF:-300}"          # Physical temperature — ALL replicas run here (K)
T_EFF_MAX="${T_EFF_MAX:-450}"  # Max effective temperature — controls λ range (K)
CUTOFF_NM="${CUTOFF_NM:-1.0}"
GAMMA_LN="${GAMMA_LN:-2.0}"
TAU_T=$(python3 -c "print(f'{1.0/float(\"${GAMMA_LN}\"):.4f}')")

# REST2 exchange & production
REPLEX_PS="${REPLEX_PS:-0.5}"
TOTAL_NS="${TOTAL_NS:-20}"
REPLEX_STEPS=$(python3 -c "print(int(round(float('${REPLEX_PS}')/float('${DT_PS}'))))")
PROD_NSTEPS=$(python3 -c "print(int(round(float('${TOTAL_NS}')*1000/float('${DT_PS}'))))")
NSTDHDL="${NSTDHDL:-500}"      # dH/dλ output interval (steps); needed for BAR/MBAR analysis

# Equilibration
EQUI_NS="${EQUI_NS:-0.2}"
EQUI_NSTEPS=$(python3 -c "print(int(round(float('${EQUI_NS}')*1000/float('${DT_PS}'))))")
NPT_SEG_STEPS="${NPT_SEG_STEPS:-10000}"
NPT_MIN_SEG="${NPT_MIN_SEG:-8}"
NPT_MAX_SEG="${NPT_MAX_SEG:-20}"
NPT_TOL_REL="${NPT_TOL_REL:-0.005}"

# Scratch — large trajectory and dhdl files go here
SCRATCH_TS="$(date +%Y%m%d-%H%M%S)"
SCRATCH_DIR="${SCRATCH_DIR:-/orcd/data/keating/001/${USER}/MD/${SLURM_JOB_ID:-$$}_${SCRATCH_TS}}"
mkdir -p "$SCRATCH_DIR"

cleanup_scratch_on_error() {
  local rc=$?
  if [[ -n "$SCRATCH_DIR" && -d "$SCRATCH_DIR" ]]; then
    if (( CURRENT_STEP >= PRESERVE_FROM_STEP )); then
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); PRESERVING scratch: $SCRATCH_DIR"
    else
      echo "[ERROR] Failed at step $CURRENT_STEP (rc=$rc); cleaning scratch: $SCRATCH_DIR"
      rm -rf "$SCRATCH_DIR"
    fi
  fi
  exit $rc
}
trap 'cleanup_scratch_on_error' ERR

SCRATCH_AVAIL=$(df -BG "$SCRATCH_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if (( SCRATCH_AVAIL < 100 )); then
  echo "[ERROR] Low disk space in scratch: ${SCRATCH_AVAIL}GB available (need ~100GB)"
  exit 1
fi

mkdir -p "${OUTDIR}"/{build,em,npt,rest2,nvt,prod,logs,trajectories} ./logs

[[ -f "$PDB_IN" ]] || { echo "[ERROR] PDB not found: $PDB_IN"; exit 1; }
echo "[OK] Parameters set  (${REPLICAS} replicas, T_REF=${T_REF} K, T_eff_max=${T_EFF_MAX} K, ${TOTAL_NS} ns)"

###############################################################################
# STEP 2 — Lambda ladder
###############################################################################
CURRENT_STEP=2
echo "===== STEP 2: REST2 lambda ladder ====="

# Effective temperatures are exponentially spaced from T_REF to T_EFF_MAX
# (same spacing as T-REMD, but only the Hamiltonian changes — not physical T).
#
#   λ_vdw_i  = 1 − T_REF/T_eff_i   ∈ [0,  1 − T_REF/T_EFF_MAX]
#   λ_coul_i = 1 − √(T_REF/T_eff_i) ∈ [0,  1 − √(T_REF/T_EFF_MAX)]

LAMBDA_ARRAYS=$(python3 - <<'PY'
import os, sys, math

n      = int(os.environ.get("REPLICAS", 48))
t_ref  = float(os.environ.get("T_REF", 300))
t_max  = float(os.environ.get("T_EFF_MAX", 450))

if n < 2:
    t_effs = [t_ref]
else:
    ratio  = (t_max / t_ref) ** (1.0 / (n - 1))
    t_effs = [t_ref * ratio**i for i in range(n)]

lam_vdw  = [max(0.0, 1.0 - t_ref / t) for t in t_effs]
lam_coul = [max(0.0, 1.0 - math.sqrt(t_ref / t)) for t in t_effs]

print("TEFF="  + ",".join(f"{t:.4f}"  for t in t_effs))
print("LVDW="  + ",".join(f"{v:.6f}"  for v in lam_vdw))
print("LCOUL=" + ",".join(f"{v:.6f}"  for v in lam_coul))
PY
)

TEFF_STR=$(echo "$LAMBDA_ARRAYS" | grep '^TEFF='  | cut -d= -f2)
VDW_STR=$( echo "$LAMBDA_ARRAYS" | grep '^LVDW='  | cut -d= -f2)
COUL_STR=$(echo "$LAMBDA_ARRAYS" | grep '^LCOUL=' | cut -d= -f2)

IFS=, read -r -a TEFF_ARR  <<< "$TEFF_STR"
IFS=, read -r -a VDW_ARR   <<< "$VDW_STR"
IFS=, read -r -a COUL_ARR  <<< "$COUL_STR"

# Space-separated strings for embedding in mdp files
VDW_MDP=$(echo  "$VDW_STR"  | tr ',' ' ')
COUL_MDP=$(echo "$COUL_STR" | tr ',' ' ')

echo "[OK] ${#VDW_ARR[@]} lambda states:"
echo "     T_eff range:   ${TEFF_ARR[0]} K  →  ${TEFF_ARR[-1]} K"
echo "     vdw-lambdas:   ${VDW_ARR[0]}  →  ${VDW_ARR[-1]}"
echo "     coul-lambdas:  ${COUL_ARR[0]} →  ${COUL_ARR[-1]}"

###############################################################################
# STEP 3 — Build system (pdb2gmx → editconf → solvate → genion)
###############################################################################
CURRENT_STEP=3
echo "===== STEP 3: Building system ====="

BUILD="${OUTDIR}/build"

echo "[INFO] Running pdb2gmx (ff=${FF}, water=${WATER})…"
$GMX pdb2gmx \
  -f "$PDB_IN" \
  -o "${BUILD}/${OUTBASE}_processed.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  -i "${BUILD}/${OUTBASE}_posre.itp" \
  -ff "$FF" -water "$WATER" -ignh \
  > "${BUILD}/pdb2gmx.log" 2>&1

echo "[INFO] Setting up ${BOX_SHAPE} box with ${BOX_BUFFER} nm buffer…"
$GMX editconf \
  -f "${BUILD}/${OUTBASE}_processed.gro" \
  -o "${BUILD}/${OUTBASE}_box.gro" \
  -bt "$BOX_SHAPE" -d "$BOX_BUFFER" -c \
  > "${BUILD}/editconf.log" 2>&1

echo "[INFO] Solvating…"
$GMX solvate \
  -cp "${BUILD}/${OUTBASE}_box.gro" \
  -cs spc216.gro \
  -o "${BUILD}/${OUTBASE}_solv.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  > "${BUILD}/solvate.log" 2>&1

cat > "${BUILD}/ions.mdp" <<EOF
integrator = steep
nsteps     = 0
EOF
$GMX grompp \
  -f "${BUILD}/ions.mdp" \
  -c "${BUILD}/${OUTBASE}_solv.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  -o "${BUILD}/ions.tpr" \
  -maxwarn 5 \
  > "${BUILD}/grompp_ions.log" 2>&1

GENION_FLAGS="-pname NA -nname CL"
(( NEUTRALIZE )) && GENION_FLAGS+=" -neutral"
python3 -c "import sys; sys.exit(0 if float('${SALT_MOLAR}')>0 else 1)" && \
  GENION_FLAGS+=" -conc ${SALT_MOLAR}"

echo "[INFO] Adding ions (neutralize=${NEUTRALIZE}, salt=${SALT_MOLAR} M)…"
echo SOL | $GMX genion \
  -s "${BUILD}/ions.tpr" \
  -o "${BUILD}/${OUTBASE}_ions.gro" \
  -p "${BUILD}/${OUTBASE}.top" \
  ${GENION_FLAGS} \
  > "${BUILD}/genion.log" 2>&1

SYS_GRO="${BUILD}/${OUTBASE}_ions.gro"
SYS_TOP="${BUILD}/${OUTBASE}.top"
[[ -f "$SYS_GRO" && -f "$SYS_TOP" ]] || { echo "[ERROR] System build failed"; exit 1; }
echo "[OK] System built"

###############################################################################
# STEP 4 — Energy minimization
###############################################################################
CURRENT_STEP=4
echo "===== STEP 4: Energy minimization ====="

EM_DIR="${OUTDIR}/em"
cat > "${EM_DIR}/em.mdp" <<EOF
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 3000

nstlist       = 10
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = ${CUTOFF_NM}
rvdw          = ${CUTOFF_NM}
pbc           = xyz

define = -DPOSRES
EOF

$GMX grompp \
  -f "${EM_DIR}/em.mdp" \
  -c "$SYS_GRO" \
  -r "$SYS_GRO" \
  -p "$SYS_TOP" \
  -o "${EM_DIR}/em.tpr" \
  -maxwarn 5 \
  > "${EM_DIR}/grompp_em.log" 2>&1

$GMX mdrun \
  -v -deffnm "${EM_DIR}/em" \
  -ntomp 8 \
  > "${EM_DIR}/mdrun_em.log" 2>&1

[[ -f "${EM_DIR}/em.gro" ]] || { echo "[ERROR] Minimization failed"; exit 1; }
echo "[OK] Minimization complete"

###############################################################################
# STEP 5 — NPT density equilibration (all at T_REF; identical logic to T-REMD)
###############################################################################
CURRENT_STEP=5
echo "===== STEP 5: NPT density equilibration (iterative) ====="

NPT_DIR="${OUTDIR}/npt"
LAST_GRO="${EM_DIR}/em.gro"
prev_vol=""
seg=0

while (( seg < NPT_MAX_SEG )); do
  seg=$((seg+1))
  echo "[INFO] NPT segment $seg"

  if (( seg == 1 )); then
    CONT_FLAGS="gen-vel = yes
gen-temp = ${T_REF}"
    SEG_DT="0.001"
  else
    CONT_FLAGS="gen-vel = no
continuation = yes"
    SEG_DT="${DT_PS}"
  fi

  cat > "${NPT_DIR}/npt_seg${seg}.mdp" <<EOF
integrator  = md
dt          = ${SEG_DT}
nsteps      = ${NPT_SEG_STEPS}
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = ${TAU_T} ${TAU_T}
ref-t       = ${T_REF} ${T_REF}

pcoupl          = C-rescale
pcoupltype      = isotropic
tau-p           = 1.0
ref-p           = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

constraints          = h-bonds
constraint-algorithm = LINCS

${CONT_FLAGS}

define = -DPOSRES
EOF

  $GMX grompp \
    -f "${NPT_DIR}/npt_seg${seg}.mdp" \
    -c "$LAST_GRO" \
    -r "$LAST_GRO" \
    -p "$SYS_TOP" \
    -o "${NPT_DIR}/npt_seg${seg}.tpr" \
    -maxwarn 5 \
    > "${NPT_DIR}/grompp_npt${seg}.log" 2>&1

  $GMX mdrun \
    -deffnm "${NPT_DIR}/npt_seg${seg}" \
    -x "${SCRATCH_DIR}/${OUTBASE}_npt_seg${seg}.xtc" \
    -ntomp 8 -nb gpu -pme gpu -bonded gpu \
    > "${NPT_DIR}/mdrun_npt${seg}.log" 2>&1

  echo Volume | $GMX energy \
    -f "${NPT_DIR}/npt_seg${seg}.edr" \
    -o "${NPT_DIR}/volume_seg${seg}.xvg" \
    > /dev/null 2>&1

  vol_avg=$(awk '!/^[#@]/ {sum+=$2; n++} END {if(n>0) printf "%.6f", sum/n}' \
            "${NPT_DIR}/volume_seg${seg}.xvg")
  echo "[INFO] Segment $seg average volume: ${vol_avg} nm^3"

  if [[ -n "$prev_vol" ]]; then
    rel_change=$(python3 -c "
prev=$prev_vol; curr=$vol_avg
print(abs(curr-prev)/prev if prev!=0 else 1.0)
")
    echo "[INFO] Relative volume change: $rel_change"
    if (( seg >= NPT_MIN_SEG )); then
      converged=$(python3 -c "print(0 if float('$rel_change')>float('${NPT_TOL_REL}') else 1)")
      if (( converged )); then
        echo "[INFO] Density converged at segment $seg"
        break
      fi
    fi
  fi

  prev_vol="$vol_avg"
  LAST_GRO="${NPT_DIR}/npt_seg${seg}.gro"
done

(( seg == NPT_MAX_SEG )) && \
  echo "[WARNING] Reached NPT_MAX_SEG=${NPT_MAX_SEG} without convergence; using last segment"

NPT_GRO="${NPT_DIR}/npt_seg${seg}.gro"
[[ -f "$NPT_GRO" ]] || { echo "[ERROR] NPT equilibration failed"; exit 1; }
echo "[OK] NPT density equilibration done"

###############################################################################
# STEP 6 — Prepare REST2 topology
###############################################################################
CURRENT_STEP=6
echo "===== STEP 6: Preparing REST2 topology ====="

REST2_DIR="${OUTDIR}/rest2"

# ── 6a: Expand all #include directives into a single flat topology ───────────
echo "[INFO] Generating preprocessed topology (expanding #includes)…"
$GMX grompp \
  -f "${EM_DIR}/em.mdp" \
  -c "$NPT_GRO" \
  -r "$NPT_GRO" \
  -p "$SYS_TOP" \
  -pp "${REST2_DIR}/${OUTBASE}_preprocessed.top" \
  -o "${REST2_DIR}/preprocess_dummy.tpr" \
  -maxwarn 10 \
  > "${REST2_DIR}/grompp_preprocess.log" 2>&1
rm -f "${REST2_DIR}/preprocess_dummy.tpr"

# ── 6b: Write topology modification script ───────────────────────────────────
cat > "${REST2_DIR}/make_rest2_top.py" <<'PYEOF'
#!/usr/bin/env python3
"""
make_rest2_top.py — Add REST2 B-state parameters to a preprocessed GROMACS topology.

For each protein atom (first non-solvent/ion moleculetype):
  • Appends a corresponding <type>_B entry to [ atomtypes ] with epsilon=0, sigma unchanged
  • Adds B-state columns to the protein's [ atoms ] section: type_B, charge_B=0, mass_B=mass_A

This enables the GROMACS free-energy (-hrex) Hamiltonian REMD exchange:
  vdw-lambdas  = 1 − T_REF/T_eff_i     → ε_protein(λ) = (1−λ)·ε_A  (scales as λ_REST2)
  coul-lambdas = 1 − √(T_REF/T_eff_i)  → q_protein(λ) = (1−λ)·q_A  (q² scales as λ_REST2)

The protein–water LJ cross-term is NOT explicitly modified; the geometric
Lorentz–Berthelot mixing rule gives sqrt(ε_protein(λ)·ε_water) = √λ_REST2·ε_cross_A
automatically.

Usage: python3 make_rest2_top.py  input_processed.top  output_rest2.top
"""
import sys, re

SOLVENT_ION_NAMES = {
    'SOL', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'SPC', 'SPCE',
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'LI', 'RB', 'CS',
    'NA+', 'CL-', 'K+', 'MG2+', 'CA2+', 'ZN2+',
}


def find_section(lines, name):
    """Return (start, end) indices for the FIRST occurrence of [ name ]."""
    target = f'[ {name} ]'
    start = None
    for i, line in enumerate(lines):
        if line.strip() == target:
            start = i
            break
    if start is None:
        return None, None
    for i in range(start + 1, len(lines)):
        if re.match(r'^\s*\[', lines[i]):
            return start, i
    return start, len(lines)


def parse_atomtypes(lines):
    """Return dict {name: sigma} from the [ atomtypes ] section."""
    s, e = find_section(lines, 'atomtypes')
    if s is None:
        return {}
    result = {}
    for line in lines[s + 1:e]:
        stripped = line.strip()
        if not stripped or stripped.startswith(';') or stripped.startswith('#'):
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        name = parts[0]
        try:
            int(parts[1])          # atomic-number column present
            sigma = float(parts[5])
        except (ValueError, IndexError):
            sigma = float(parts[4])
        result[name] = sigma
    return result


def get_molecule_ranges(lines):
    """Return [(mol_name, mol_start, mol_end), ...] for all moleculetypes."""
    system_start = next(
        (i for i, l in enumerate(lines)
         if l.strip() in ('[ system ]', '[ molecules ]')),
        len(lines),
    )
    mol_starts = [i for i, l in enumerate(lines) if l.strip() == '[ moleculetype ]']
    ranges = []
    for k, ms in enumerate(mol_starts):
        me = mol_starts[k + 1] if k + 1 < len(mol_starts) else system_start
        name = None
        for j in range(ms + 1, me):
            s = lines[j].strip()
            if s and not s.startswith(';') and not s.startswith('[') and not s.startswith('#'):
                name = s.split()[0]
                break
        ranges.append((name, ms, me))
    return ranges


def modify_topology(input_path, output_path):
    with open(input_path) as f:
        lines = f.readlines()

    # ── Identify protein molecule ────────────────────────────────────────────
    protein_range = None
    for name, ms, me in get_molecule_ranges(lines):
        if name and name.upper() not in SOLVENT_ION_NAMES:
            protein_range = (name, ms, me)
            break
    if protein_range is None:
        sys.exit('[ERROR] make_rest2_top: no protein molecule found in topology')

    prot_name, prot_s, prot_e = protein_range
    print(f'[REST2] Protein molecule: "{prot_name}" (lines {prot_s}–{prot_e})',
          file=sys.stderr)

    # ── Find [ atoms ] within protein molecule ───────────────────────────────
    prot_atoms_s = None
    for i in range(prot_s, prot_e):
        if lines[i].strip() == '[ atoms ]':
            prot_atoms_s = i
            break
    if prot_atoms_s is None:
        sys.exit('[ERROR] make_rest2_top: no [ atoms ] section in protein molecule')

    prot_atoms_e = prot_e
    for i in range(prot_atoms_s + 1, prot_e):
        if re.match(r'^\s*\[', lines[i]):
            prot_atoms_e = i
            break

    # ── Collect unique protein atom types ────────────────────────────────────
    protein_types = set()
    for line in lines[prot_atoms_s + 1:prot_atoms_e]:
        s = line.strip()
        if s and not s.startswith(';') and not s.startswith('#'):
            parts = s.split()
            if len(parts) >= 7:
                protein_types.add(parts[1])

    print(f'[REST2] {len(protein_types)} unique protein atom types: '
          f'{sorted(protein_types)}', file=sys.stderr)

    # ── Parse sigma values from [ atomtypes ] ────────────────────────────────
    atomtypes_sigma = parse_atomtypes(lines)
    atomtypes_s, atomtypes_e = find_section(lines, 'atomtypes')

    # ── Build modified topology (two-pass write) ─────────────────────────────
    out = []
    i = 0
    while i < len(lines):
        # After [ atomtypes ] block: inject _B type definitions
        if atomtypes_e is not None and i == atomtypes_e:
            out.append('; -------- REST2 B-state atom types (epsilon = 0) --------\n')
            for ptype in sorted(protein_types):
                sigma = atomtypes_sigma.get(ptype, 3.5e-1)
                out.append(
                    f'{ptype + "_B":<14}  0.00000  0.00000  A'
                    f'  {sigma:14.8e}  {0.0:14.8e}\n'
                )
            out.append('; ---------------------------------------------------------\n\n')
            atomtypes_e = None  # only inject once

        # Within protein [ atoms ]: add B-state columns
        if i == prot_atoms_s:
            out.append(lines[i])   # '[ atoms ]' header
            i += 1
            while i < prot_atoms_e:
                aline = lines[i]
                s = aline.strip()
                if s and not s.startswith(';') and not s.startswith('#'):
                    parts = s.split()
                    if len(parts) >= 7:
                        # Standard columns: nr type resnr resname atomname cgnr charge [mass]
                        nr, atype, resnr, resname, atomname, cgnr, charge = parts[:7]
                        mass = parts[7] if len(parts) >= 8 else '0.0'
                        b_type = atype + '_B'
                        out.append(
                            f'{nr:>6} {atype:>10} {resnr:>5} {resname:>6} {atomname:>6}'
                            f' {cgnr:>5} {charge:>12} {mass:>12}'
                            f'  {b_type:>12} {"0.000000":>12} {mass:>12}\n'
                        )
                        i += 1
                        continue
                out.append(aline)
                i += 1
            continue   # i already advanced inside inner loop

        out.append(lines[i])
        i += 1

    with open(output_path, 'w') as f:
        f.writelines(out)
    print(f'[REST2] Modified topology written: {output_path}', file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f'Usage: {sys.argv[0]} input_processed.top output_rest2.top')
    modify_topology(sys.argv[1], sys.argv[2])
PYEOF

# ── 6c: Run topology modification ────────────────────────────────────────────
echo "[INFO] Adding REST2 B-state parameters to protein atoms…"
python3 "${REST2_DIR}/make_rest2_top.py" \
  "${REST2_DIR}/${OUTBASE}_preprocessed.top" \
  "${REST2_DIR}/${OUTBASE}_rest2.top"

# All downstream steps use this topology
REST2_TOP="${REST2_DIR}/${OUTBASE}_rest2.top"
[[ -f "$REST2_TOP" ]] || { echo "[ERROR] REST2 topology generation failed"; exit 1; }
echo "[OK] REST2 topology ready: $REST2_TOP"

###############################################################################
# STEP 7 — Prepare per-replica NVT equilibration directories
###############################################################################
CURRENT_STEP=7
echo "===== STEP 7: Preparing per-replica NVT equilibration ====="

NVT_DIR="${OUTDIR}/nvt"

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  REPDIR="${NVT_DIR}/rep${idx}"
  mkdir -p "$REPDIR"

  cat > "${REPDIR}/nvt.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${EQUI_NSTEPS}
nstxout-compressed = 2500
nstlog      = 1000
nstenergy   = 1000

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl  = V-rescale
tc-grps = Protein Non-Protein
tau-t   = ${TAU_T} ${TAU_T}
ref-t   = ${T_REF} ${T_REF}

pcoupl = no

constraints          = h-bonds
constraint-algorithm = LINCS

gen-vel  = yes
gen-temp = ${T_REF}

nstcomm = 500

; ── REST2 free-energy settings ─────────────────────────────────────────────
; Each replica equilibrates in its own scaled Hamiltonian.
; No exchange during NVT — only the Hamiltonian/lambda state varies.
free-energy           = yes
init-lambda-state     = ${i}
vdw-lambdas           = ${VDW_MDP}
coul-lambdas          = ${COUL_MDP}
sc-alpha              = 0            ; linear mixing (no soft-core)
calc-lambda-neighbors = 1            ; compute ±1 neighbor energies
nstdhdl               = ${NSTDHDL}
separate-dhdl-file    = yes
EOF

  $GMX grompp \
    -f "${REPDIR}/nvt.mdp" \
    -c "$NPT_GRO" \
    -r "$NPT_GRO" \
    -p "$REST2_TOP" \
    -o "${REPDIR}/nvt.tpr" \
    -maxwarn 10 \
    > "${REPDIR}/grompp.log" 2>&1

  # Pre-create symlink so deffnm writes xtc directly to scratch
  ln -sf "${SCRATCH_DIR}/nvt_rep${idx}.xtc" "${REPDIR}/nvt.xtc"
done
echo "[OK] NVT .tpr files prepared for all replicas"

###############################################################################
# STEP 8 — Run per-replica NVT equilibration
###############################################################################
CURRENT_STEP=8
echo "===== STEP 8: Running NVT equilibration ====="

MULTIDIR_NVT=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  MULTIDIR_NVT+="${NVT_DIR}/rep${idx} "
done

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir $MULTIDIR_NVT \
  -deffnm nvt \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_nvt_equil.log" 2>&1

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  [[ -f "${NVT_DIR}/rep${idx}/nvt.gro" ]] || \
    { echo "[ERROR] NVT equil failed for replica $idx"; exit 1; }
done
echo "[OK] NVT equilibration done"

###############################################################################
# STEP 9 — Prepare REST2 production inputs
###############################################################################
CURRENT_STEP=9
echo "===== STEP 9: Preparing REST2 production inputs ====="

PROD_DIR="${OUTDIR}/prod"
NSTXOUT_COMP=$(python3 -c "print(int(${REPLEX_STEPS}*10))")

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  REPDIR="${PROD_DIR}/rep${idx}"
  mkdir -p "$REPDIR"

  cat > "${REPDIR}/remd.mdp" <<EOF
integrator  = md
dt          = ${DT_PS}
nsteps      = ${PROD_NSTEPS}
nstxout-compressed = ${NSTXOUT_COMP}
nstlog      = ${NSTXOUT_COMP}
nstenergy   = ${NSTXOUT_COMP}

cutoff-scheme = Verlet
nstlist     = 10
coulombtype = PME
rcoulomb    = ${CUTOFF_NM}
rvdw        = ${CUTOFF_NM}
DispCorr    = EnerPres
pbc         = xyz

tcoupl  = V-rescale
tc-grps = Protein Non-Protein
tau-t   = ${TAU_T} ${TAU_T}
ref-t   = ${T_REF} ${T_REF}

pcoupl = no

constraints          = h-bonds
constraint-algorithm = LINCS

gen-vel      = no
continuation = yes

nstcomm = 500

; ── REST2 free-energy / Hamiltonian REMD settings ──────────────────────────
; vdw-lambdas / coul-lambdas are the same for all replicas.
; Only init-lambda-state differs — it selects which λ value this replica runs at.
; calc-lambda-neighbors = 1 lets each replica evaluate energy at adjacent λ
; states, which is what -hrex needs to compute the exchange criterion.
free-energy           = yes
init-lambda-state     = ${i}
vdw-lambdas           = ${VDW_MDP}
coul-lambdas          = ${COUL_MDP}
sc-alpha              = 0
calc-lambda-neighbors = 1
nstdhdl               = ${NSTXOUT_COMP}
separate-dhdl-file    = yes
EOF

  $GMX grompp \
    -f "${REPDIR}/remd.mdp" \
    -c "${NVT_DIR}/rep${idx}/nvt.gro" \
    -t "${NVT_DIR}/rep${idx}/nvt.cpt" \
    -p "$REST2_TOP" \
    -o "${REPDIR}/remd.tpr" \
    -maxwarn 10 \
    > "${REPDIR}/grompp.log" 2>&1

  # Pre-create symlinks for large files on scratch
  ln -sf "${SCRATCH_DIR}/remd_rep${idx}.xtc"  "${REPDIR}/remd.xtc"
  ln -sf "${SCRATCH_DIR}/remd_rep${idx}.dhdl.xvg" "${REPDIR}/remd.dhdl.xvg"
done
echo "[OK] REST2 production .tpr files prepared"

###############################################################################
# STEP 10 — Run REST2 production
###############################################################################
CURRENT_STEP=10
echo "===== STEP 10: Running REST2 production ====="
echo "[INFO] ${REPLICAS} replicas × ${TOTAL_NS} ns, exchange every ${REPLEX_PS} ps"
echo "[INFO] All replicas at T=${T_REF} K; effective T range ${TEFF_ARR[0]}–${TEFF_ARR[-1]} K"

MULTIDIR_PROD=""
for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  MULTIDIR_PROD+="${PROD_DIR}/rep${idx} "
done

mpirun -np "$REPLICAS" "$GMX" mdrun \
  -multidir $MULTIDIR_PROD \
  -deffnm remd \
  -replex "$REPLEX_STEPS" \
  -ntomp 1 \
  -nb gpu -pme gpu -bonded gpu \
  -pin on \
  > "${OUTDIR}/logs/mdrun_rest2.log" 2>&1

REMD_EXIT=$?
if (( REMD_EXIT != 0 )); then
  echo "[ERROR] REST2 production failed (exit=$REMD_EXIT)"
  exit 1
fi
echo "[OK] REST2 production done"

###############################################################################
# STEP 11 — Finalize outputs & symlinks
###############################################################################
CURRENT_STEP=11
echo "===== STEP 11: Finalizing outputs ====="

for i in $(seq 0 $((REPLICAS-1))); do
  idx=$(printf "%03d" "$i")
  for prefix in nvt remd; do
    src="${SCRATCH_DIR}/${prefix}_rep${idx}.xtc"
    [[ -f "$src" ]] && ln -sf "$src" "${OUTDIR}/trajectories/${prefix}_rep${idx}.xtc"
  done
  # dhdl files (used for BAR/MBAR free-energy analysis across lambda states)
  src="${SCRATCH_DIR}/remd_rep${idx}.dhdl.xvg"
  [[ -f "$src" ]] && ln -sf "$src" "${OUTDIR}/trajectories/remd_rep${idx}.dhdl.xvg"
done

for f in "${SCRATCH_DIR}/${OUTBASE}_npt_seg"*.xtc; do
  [[ -f "$f" ]] && ln -sf "$f" "${OUTDIR}/trajectories/$(basename "$f")"
done

echo 0 | $GMX trjconv \
  -s "${PROD_DIR}/rep000/remd.tpr" \
  -f "${PROD_DIR}/rep000/remd.gro" \
  -o "${OUTDIR}/${OUTBASE}_final_rep000.pdb" \
  -pbc whole \
  > /dev/null 2>&1 || true

echo "[OK] Outputs finalized"

###############################################################################
# STEP 12 — Write parameters log
###############################################################################
CURRENT_STEP=12
ELAPSED_S=$(( $(date +%s) - CURRENT_TIME ))

for dest in "${OUTDIR}/parameters.txt" "${SCRATCH_DIR}/parameters.txt"; do
  cat > "$dest" <<EOF
Elapsed seconds:         $ELAPSED_S
Job Name:                ${SLURM_JOB_NAME:-rest2}
Job ID:                  ${SLURM_JOB_ID:-N/A}
Method:                  REST2 (Replica Exchange with Solute Tempering 2)
Replicas:                $REPLICAS
Total ns per replica:    $TOTAL_NS
Physical temperature:    ${T_REF} K  (all replicas)
Max effective temp:      ${T_EFF_MAX} K
Lambda vdw range:        ${VDW_ARR[0]}  →  ${VDW_ARR[-1]}
Lambda coul range:       ${COUL_ARR[0]} →  ${COUL_ARR[-1]}
Force field:             $FF
Water model:             $WATER
Box shape:               $BOX_SHAPE
Box buffer (nm):         $BOX_BUFFER
Neutralize:              $NEUTRALIZE
Salt molarity:           $SALT_MOLAR
Cutoff (nm):             $CUTOFF_NM
tau_t (ps):              $TAU_T
Time step (ps):          $DT_PS
Exchange interval (ps):  $REPLEX_PS
Exchange interval (steps): $REPLEX_STEPS
Production steps:        $PROD_NSTEPS
dH/dλ output (steps):   $NSTDHDL
Output directory:        $OUTDIR
Scratch directory:       $SCRATCH_DIR
REST2 topology:          $REST2_TOP
GPUs:                    ${SLURM_GPUS_ON_NODE:-unknown}
NPT convergence segs:    $seg
EOF
done

###############################################################################
# STEP 13 — Summary
###############################################################################
echo ""
echo "=========================================="
echo "===== JOB COMPLETED $(date) ====="
echo "=========================================="
echo "[INFO] Output directory:    $OUTDIR"
echo "[INFO] Trajectory links:    ${OUTDIR}/trajectories/"
echo "[INFO] Actual trajectories: $SCRATCH_DIR"
echo "[INFO] Exchange log:        ${PROD_DIR}/rep000/remd.log"
echo ""
echo "Next steps:"
echo "  1. Check exchange acceptance:"
echo "       grep 'Repl  average probabilities' ${PROD_DIR}/rep000/remd.log"
echo ""
echo "  2. Analyze reference (λ=0) replica at T_REF:"
echo "       trajectories/remd_rep000.xtc  +  ${PROD_DIR}/rep000/remd.tpr"
echo ""
echo "  3. Demux to continuous-lambda trajectories (Hamiltonian walk):"
echo "       $GMX demux ${PROD_DIR}/rep000/remd.log"
echo ""
echo "  4. Free-energy analysis across lambda states (BAR/MBAR):"
echo "       $GMX bar -f trajectories/remd_rep*.dhdl.xvg"
echo "       or use alchemical-analysis / pymbar on the dhdl.xvg files"
echo "=========================================="
