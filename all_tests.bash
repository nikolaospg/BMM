#!/bin/bash -e
#This script is used to run the source files that we have created for test purposes. They are found in the /test folder.
#The first command line argument on the product_test is the n dimension, and the second is a flag specifying whether the matrix is done with a Filter or not.
#The command line argument on the conversion_test is the n dimension.
#The number of iterations is specified in the for statement.

#The way you are supposed to use this script is by simply running the make test command. It creates the binary files on the .bin folder and executes them
#with the help of this script.

for (( i=0 ; i<10; i++))
do
  ./bin/product_test 500 0.9 ss 0
  ./bin/product_test 500 0.9 ss 1
  ./bin/product_test 500 0.9 sp 0
  ./bin/product_test 50 0.9 bs 0 25
  ./bin/product_test 50 0.9 bs 1 25
  ./bin/product_test 20 0.9 cp 1
  ./bin/conversion_test 500
  ./bin/blocking_test 500 20
  ./bin/temp_test 5000 2635
  echo "Iteration $i Complete - No Errors"
done 
