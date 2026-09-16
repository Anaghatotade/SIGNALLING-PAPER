#!/bin/bash
#SBATCH --job-name=N_250 # Job name
#SBATCH --output=lammps_test.out        # Output file
#SBATCH --error=lammps_test.err         # Error file
#SBATCH --ntasks=64                      # Number of MPI tasks (processes)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --partition=compute

# Load necessary modules
source /apps/compilers/intel/oneapi/setvars.sh
module load lammps-cpu

# Change to the directory where the job was submitted from
cd $SLURM_SUBMIT_DIR

export I_MPI_FABRICS=shm:tcp
export I_MPI_OFI_PROVIDER=tcp


mpirun -np 64 /home/dr_gaurav_chauhan1/lammps_new/build/lmp -in lammps_sim.in
