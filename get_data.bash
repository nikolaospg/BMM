#!/bin/bash

#This script helps us download and extract the matrix market files, specified on the dataset.txt file. They get stored on the data directory.
#The way you are supposed to use this script is by simply running the make data command.
mkdir -p data
cat dataset.txt | xargs wget -P data
cd data
for file in ./*.tar*
do
    tar -xzf ${file} --strip-components 1
done
rm  ./*.tar*
