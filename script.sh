#!/bin/bash -e
#bash -e script.sh -> Recommended syntax
#Script used to make tests on whether the programms are correct. Just uncomment the test programm you want for it to run
for (( i=0 ; i<5000; i++))
do
 # ./convertion_test.out 8000
  ./testproduct.out 500 0
  echo "$i"
done 
