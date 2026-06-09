#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p pi_keating
#SBATCH --mem=80000
#SBATCH --nodelist=node[3619-3620]
#SBATCH -o ./logs/install_gromacs-plumed%j.out
#SBATCH -e ./logs/install_gromacs-plumed%j.err
#SBATCH --constraint="rocky8"
#SBATCH --gres=gpu:1


# module unload miniforge/23.11.0-0
# openmpi/4.1.4
export LD_LIBRARY_PATH=
# module load nvhpc/26.1
module load openmpi/5.0.8
module load cuda/12.9.1

# export PLUMED_KERNEL=/home/jhalpin/opt/plumed/2.9.4/lib/libplumedKernel.so
# export LD_LIBRARY_PATH=/home/jhalpin/opt/plumed/2.9.4/lib:$LD_LIBRARY_PATH
source ~/plumed.sh

cd ./gromacs-2024.3/
rm -rf build
mkdir build
plumed patch -p -e gromacs-2024.3
cd build

cmake .. -DGMX_THREAD_MPI=OFF -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DGMX_MPI=on -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs/2024.3-plumed
# cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DGMX_MPI=on -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs/2024.3 -DGMX_NVSHMEM=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
# cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DGMX_MPI=on -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs/2024.2
make -j 8
make check
make install

# #SBATCH --nodelist=node3500
