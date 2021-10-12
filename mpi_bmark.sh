#!/bin/bash
#SBATCH --time=00:10:00               # Run time (days-hh:mm:ss) - (max 7days) 
#SBATCH --partition=batch             # Submit to queue/partition named batch
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load gcc/9.2.0 openmpi/3.1.4

make mpi
make bash
echo "Starting benchmark MPI"
./bmm_self_mpi.bash data/europe_osm.mtx 4
