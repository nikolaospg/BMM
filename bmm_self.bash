#!/bin/bash
# Utility script to run BMM with A=B(=F)
# Usage: ./bmm_self <matrix market file> <method> [-f]

filt=""
if [ "$3" = "-f" ]; then
    filt="-f $1"
fi

./bin/bmm -a $1 -b $1 $filt -m $2
