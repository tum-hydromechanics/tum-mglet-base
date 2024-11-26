#!/bin/bash
#SBATCH -J MGLET
#SBATCH --account=pn52gi
#SBATCH --time=0:01:00
#SBATCH --export=NONE
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#SBATCH -D ./

export OMP_NUM_THREADS=1

module load slurm_setup
module load hdf5/1.14.3-intel24-impi
module switch mpi_settings/2.0 mode=cpu-only

srun -n 4 ./mglet 
