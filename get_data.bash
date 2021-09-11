#!/bin/bash
mkdir -p data
cat dataset.txt | xargs wget -P data
cd data
for file in ./*.tar*
do
    tar -xzf ${file} --strip-components 1
done
rm  ./*.tar*
