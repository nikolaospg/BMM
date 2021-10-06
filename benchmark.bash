#!/bin/bash
# Run BMM on some matrices with a given method
# Usage: ./benchmark <method> [-f]

for f in belgium_osm.mtx com-Youtube.mtx dblp-2010.mtx delaunay_n24.mtx europe_osm.mtx hollywood-2009.mtx mycielskian13.mtx
do
    ./bmm_self.bash data/$f $1 $2
done
