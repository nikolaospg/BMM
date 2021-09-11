#!/bin/bash -e
#bash -e script.sh -> Recommended syntax

for (( i=0 ; i<10; i++))
do
  ./bin/product_test 500 0
  ./bin/product_test 500 1
  ./bin/conversion_test 500
  echo "Iteration $i"
done 
