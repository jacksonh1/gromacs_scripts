#!/bin/bash
# VALIDATION ONLY — check that REST2 topology scaling produces the right ENERGIES for a
# force field, before trusting REST2 runs with it. CPU-only, a few minutes; no GPU needed.
#
#   bash scripts/validation/validate_rest2_scaling.sh SYSTEM.gro SYSTEM.top WORKDIR [LAMBDA] [TRAJ]
#
# SYSTEM.gro/.top: a solvated, ionized system as the REST2 engine builds it (build/ or em/
# of a job, or pdb2gmx+solvate+genion by hand). TRAJ (optional .xtc/.trr/.gro) supplies
# several frames; default is SYSTEM.gro alone. LAMBDA defaults to 0.5. Three checks, each a
# frame-by-frame per-term ratio against the UNSCALED topology (rest2_energy_ratios.py):
#
#   (a) protein marked, lambda=1.0  -> every term x1 (partial_tempering at 1.0 is a no-op)
#   (b) ALL atoms marked, lambda    -> bonds/angles x1, every other term xlambda
#   (c) protein marked, lambda      -> [energy groups] protein-protein xlambda,
#                                      protein-water x sqrt(lambda), water-water x1,
#                                      torsions and 1-4 xlambda, bonds/angles x1
#
# plus check_rest2_topology.py on the (c) pair. Exit 1 if any check fails.
# See docs/FORCE_FIELDS.md ("Validating a new force field for REST2").
set -euo pipefail

GRO="$(realpath "$1")"; TOP="$(realpath "$2")"; WORK="$3"; LAM="${4:-0.5}"; TRAJ="${5:-$GRO}"
TRAJ="$(realpath "$TRAJ")"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SIM="${REPO}/scripts/simulation"
VAL="${REPO}/scripts/validation"

source "${REPO}/site_config.sh"
# same environment as REST2-gromacs.sbatch STEP 0
module purge
module load gcc
module load deprecated-modules          # provides REST2_CUDA_MODULE (cuda 12.4)
module load "$REST2_CUDA_MODULE"
module load "$OPENMPI_MODULE"
set +u; source "$REST2_GMXRC"; source "$PLUMED_SH"; set -u
export LD_LIBRARY_PATH="${HCOLL_COMPAT_DIR}:${LD_LIBRARY_PATH:-}"
GMX="${GMX:-gmx_mpi}"

mkdir -p "$WORK"; cd "$WORK"
# grompp resolves the .top's relative #includes (posre.itp, *_chain.itp) from its own dir
TOPDIR="$(dirname "$TOP")"

cat > pp.mdp <<EOF
integrator    = md
nsteps        = 0
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
EOF
cat > rerun.mdp <<EOF
integrator    = md
nsteps        = 0
continuation  = yes
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
DispCorr      = EnerPres
energygrps    = Protein non-Protein
EOF

echo "[INFO] processed topology (no -DPOSRES, as STEP 7)…"
( cd "$TOPDIR" && $GMX grompp -f "$WORK/pp.mdp" -c "$GRO" -p "$TOP" -pp "$WORK/processed.top" \
    -o "$WORK/pp.tpr" -maxwarn 5 ) > grompp_pp.log 2>&1
python3 "${SIM}/mark_hot_region.py" processed.top > marked.top
python3 "${VAL}/mark_all_atoms.py" processed.top > allhot.top
plumed partial_tempering 1.0    < marked.top > prot_1.top    2> pt_prot_1.log
plumed partial_tempering "$LAM" < marked.top > prot_lam.top  2> pt_prot_lam.log
plumed partial_tempering "$LAM" < allhot.top > all_lam.top   2> pt_all_lam.log

python3 "${SIM}/check_rest2_topology.py" prot_1.top prot_lam.top --lambda-cold 1.0 --lambda-hot "$LAM"

rerun() {  # rerun <top> <name>: single-point energies over TRAJ -> <name>.xvg (all terms but Potential)
  $GMX grompp -f rerun.mdp -c "$GRO" -p "$1" -o "$2.tpr" -maxwarn 2 > "grompp_$2.log" 2>&1
  $GMX mdrun -s "$2.tpr" -deffnm "$2" -rerun "$TRAJ" -nb cpu -ntomp "${NTOMP:-8}" > "mdrun_$2.log" 2>&1
  local listing sel
  # `gmx energy` prints its term table, then aborts on the empty selection "0": the non-zero
  # exit is expected here, the table is all we want.
  listing="$(echo 0 | $GMX energy -f "$2.edr" 2>&1 || true)"
  sel="$(grep -E '^( +[0-9]+ +[^ ]+)+ *$' <<< "$listing" | grep -oE '[0-9]+ +[^ ]+' \
         | awk '$2 != "Potential" {print $1}' | tr '\n' ' ')"
  [[ -n "$sel" ]] || { echo "[ERROR] no energy terms listed for $2" >&2; exit 1; }
  echo "$sel 0" | $GMX energy -f "$2.edr" -o "$2.xvg" > "energy_$2.log" 2>&1
}
echo "[INFO] reruns (CPU)…"
rerun processed.top ref
rerun prot_1.top    prot_1
rerun prot_lam.top  prot_lam
rerun all_lam.top   all_lam

# Totals whose pairs are all hot or all cold get a single ratio; mixed totals are '*' and
# their energy-group parts are checked instead.
TERMS_ALL1="Bond=1,Angle=1,Proper Dih.=1,Per. Imp. Dih.=1,LJ-14=1,Coulomb-14=1,LJ (SR)=1,Disper. corr.=1,Coulomb (SR)=1,Coul. recip.=1"
GROUPS_ALL1="Coul-SR:Protein-Protein=1,LJ-SR:Protein-Protein=1,Coul-14:Protein-Protein=1,LJ-14:Protein-Protein=1,Coul-SR:Protein-non-Protein=1,LJ-SR:Protein-non-Protein=1,Coul-14:Protein-non-Protein=1,LJ-14:Protein-non-Protein=1,Coul-SR:non-Protein-non-Protein=1,LJ-SR:non-Protein-non-Protein=1,Coul-14:non-Protein-non-Protein=1,LJ-14:non-Protein-non-Protein=1"
L="$LAM"; S="sqrt($LAM)"
EXP_B="Bond=1,Angle=1,Proper Dih.=$L,Per. Imp. Dih.=$L,LJ-14=$L,Coulomb-14=$L,LJ (SR)=$L,Disper. corr.=$L,Coulomb (SR)=$L,Coul. recip.=$L,$(echo "$GROUPS_ALL1" | sed "s/=1/=$L/g")"
EXP_C="Bond=1,Angle=1,Proper Dih.=$L,Per. Imp. Dih.=$L,LJ-14=$L,Coulomb-14=$L,LJ (SR)=*,Disper. corr.=*,Coulomb (SR)=*,Coul. recip.=*,Coul-SR:Protein-Protein=$L,LJ-SR:Protein-Protein=$L,Coul-14:Protein-Protein=$L,LJ-14:Protein-Protein=$L,Coul-SR:Protein-non-Protein=$S,LJ-SR:Protein-non-Protein=$S,Coul-14:Protein-non-Protein=$S,LJ-14:Protein-non-Protein=$S,Coul-SR:non-Protein-non-Protein=1,LJ-SR:non-Protein-non-Protein=1,Coul-14:non-Protein-non-Protein=1,LJ-14:non-Protein-non-Protein=1"

rc=0
echo; echo "===== (a) protein marked, lambda=1.0: must equal the unscaled topology ====="
python3 "${VAL}/rest2_energy_ratios.py" ref.xvg prot_1.xvg "${TERMS_ALL1},${GROUPS_ALL1}" --tol 1e-6 || rc=1
echo; echo "===== (b) all atoms marked, lambda=${LAM}: bonds/angles x1, everything else x${LAM} ====="
python3 "${VAL}/rest2_energy_ratios.py" ref.xvg all_lam.xvg "$EXP_B" || rc=1
echo; echo "===== (c) protein marked, lambda=${LAM}: the REST2 split ====="
python3 "${VAL}/rest2_energy_ratios.py" ref.xvg prot_lam.xvg "$EXP_C" || rc=1
echo
if (( rc )); then echo "[FAIL] REST2 scaling validation failed — see above"; else echo "[OK] REST2 scaling validated (lambda=${LAM})"; fi
exit $rc
