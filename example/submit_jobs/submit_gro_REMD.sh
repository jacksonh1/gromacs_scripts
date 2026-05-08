
# ---------------------------------------------
REPLICAS=48
T_MAX=450
TOTAL_NS=2
REPLEX_PS=1



PDB_IN="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/input_pdbs/helix_fusion.pdb"
OUTBASE="helix_fusion"
OUTDIR="/home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/example/outputs/output_T-REMD/"$OUTBASE"-"$TOTAL_NS"ns-REMD-300-"$T_MAX"K-"$REPLICAS"reps-NVT-exf-"$REPLEX_PS"ps"

sbatch -n "$REPLICAS" --export=ALL,REPLEX_PS=$REPLEX_PS,TOTAL_NS=$TOTAL_NS,T_MAX=$T_MAX,PDB_IN=$PDB_IN,OUTDIR=$OUTDIR,OUTBASE=$OUTBASE /home/jhalpin/orcd/pool/09-fragfold/RELE_simulations/gromacs_REMD/gromacs_scripts/REMD-gromacs.sbatch

