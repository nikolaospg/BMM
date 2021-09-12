#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>


/** 
 * Prints a vector(row major fashion)
 * */
void print_vector(int* vector, int size){
    for(int i=0; i<size; i++){
        printf("%d  ", vector[i]);
    }   
    printf("\n\n");
}

// Prints a vector as 2D matrix (row-major)
void print_vector2D(int* vector, int rows, int cols){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            printf("%d  ", vector[i*cols+j]);
        }
        printf("\n");
    }   
    printf("\n");
}

/*Creates a random vector with boolean values 
 * input: The size of the vector
 * returns the pointer.*/
int* random_vector(int size){
    int* ret=(int*)malloc(size*sizeof(int*));
    if(ret==NULL){
        printf("Could not allocate memory in random_vector function. Exiting\n");
        exit(-1);
    }   

    double temp;
    for(int i=0; i<size; i++){
        temp=rand()/(double)RAND_MAX;
        if(temp>0.5){           //The way I convert this probability distribution to a discrete one, is by getting the ones >0.5 as equal to 1, and the rest as equal to 0
            ret[i]=1;
        }
        else{
            ret[i]=0;
        }
    }   
    return ret;
}

//Returns 1 if the requested element exists in the array it returns 0.
//Args: int* a-> Pointer to the array
//      initial_index ->    The index of the first element of the array (the user just passes 0), it is only different for the recursive calls inside the function.
//      size->          The number of elements of the array
//      elem->        The wanted element
//      The reason we pass the initial_index argument is to be able to know the actual index of the element we get. The algorithm works by cutting the array and passing the parts of it
//      as arguments, so in the recursive calls the information about the initial index is lost.
int binary_search(int* a,int begin, int end, int elem){
    if(begin <= end){
        int mid_idx = begin + (end-begin)/2;
        if(a[mid_idx]==elem){
            return 1;
        } else if(a[mid_idx]>elem){
            return binary_search(a, begin, mid_idx-1, elem);
        } else{
            return binary_search(a, mid_idx+1, end, elem);
        }
    }   
    
    return 0;
}
