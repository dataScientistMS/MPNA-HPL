#!/usr/bin/env bash
#SBATCH --account=r250142
#SBATCH --constraint=x64cpu
#SBATCH --partition=short
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --job-name=HPL_S1_N_100000
#SBATCH --output=job.out
#SBATCH --error=job.err

romeo_load_x64cpu_env
module load aocl-4.2.0 openmpi/aocc-4.2.0/4.1.7

export OMP_NUM_THREADS=1
export BLIS_JC_NT=1
export BLIS_IC_NT=1
export BLIS_JR_NT=1
export BLIS_IR_NT=1


cd $SLURM_SUBMIT_DIR
srun ../xhpl
