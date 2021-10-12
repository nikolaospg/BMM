#!/bin/bash
#SBATCH --time=00:30:00               # Run time (days-hh:mm:ss) - (max 7days) 
#SBATCH --partition=batch             # Submit to queue/partition named batch
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=10

make mpi
make bash

module load gcc/9.2.0 openmpi/3.1.4
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mat="data/mycielskian13.mtx"
min_nb=2
max_nb=7
nb=$min_nb

echo "Starting benchmark MPI big"
while [[ $nb -le $max_nb ]]
do
    ./bmm_self_mpi.bash $mat nb
    ./bmm_self_mpi.bash $mat nb -f
    ((nb=nb+1))
done

mat="data/belgium_osm.mtx"
min_nb=2
max_nb=7
nb=$min_nb

echo "Starting benchmark MPI big"
while [[ $nb -le $max_nb ]]
do
    ./bmm_self_mpi.bash $mat $nb
    ./bmm_self_mpi.bash $mat $nb -f
    ((nb=nb+1))
done
