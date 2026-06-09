#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p pi_keating
#SBATCH --mem=80000
#SBATCH -o ./logs/install_plumed%j.out
#SBATCH -e ./logs/install_plumed%j.err
#SBATCH --nodelist=node[3619-3620]
#SBATCH --constraint="rocky8"


module purge
export LD_LIBRARY_PATH=
module load openmpi/5.0.8
module load cuda/12.9.1
cd ./plumed-2.9.4/
make distclean 2>/dev/null || true
./configure --prefix=$HOME/opt/plumed/2.9.4
make -j 8
make regtest
make install

# #SBATCH --nodelist=node3500
