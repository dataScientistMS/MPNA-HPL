#!/bin/bash
#SBATCH --account=r250142
#SBATCH --constraint=x64cpu
#SBATCH --partition=short
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --job-name=HPL_CUSTOM
#SBATCH --output=custom_out_%j.out

# Rechargement des modules
romeo_load_x64cpu_env
module load aocl-4.2.0 openmpi/aocc-4.2.0/4.1.7

# Lancement avec arguments N=115000 NB=224
srun ./xhpl_custom 115000 224
