#!/bin/bash
# Utility script to run MPI BMM with A=B(=F)
# Usage: ./bmm_self_mpi <matrix market file> <nblocks> [-f]

filt=""
if [ "$3" = "-f" ]; then
    filt="-f $1"
fi

srun bin/mpi -a $1 -b $1 -s $2 $filt
