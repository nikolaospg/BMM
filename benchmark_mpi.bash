#!/bin/bash
# Run BMM on some matrices with a given method
# Usage: ./benchmark <nblocks> [-f]

for f in belgium_osm.mtx com-Youtube.mtx dblp-2010.mtx mycielskian13.mtx com-Youtube.mtx delaunay_n24.mtx europe_osm.mtx hollywood-2009.mtx
do
    ./bmm_self_mpi.bash data/$f $1 $2
done
