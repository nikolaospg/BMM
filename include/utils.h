#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>


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
}

/*Creates a random vector with boolean values 
 * input: The size of the vector
 * sparsity: The sparsity of the vector. Each element will have a sparsity chance of being a zero
 * returns the pointer.*/
int* random_vector(int size, double sparsity){
    int* ret=(int*)malloc(size*sizeof(int*));
    if(ret==NULL){
        printf("Could not allocate memory in random_vector function. Exiting\n");
        exit(-1);
    }
    assert(sparsity >=0 && sparsity <=1);

    double temp;
    for(int i=0; i<size; i++){
        temp=rand()/(double)RAND_MAX;
        if(temp>sparsity){
            ret[i]=1;
        } else {
            ret[i]=0;
        }
    }   

    return ret;
}

//Returns 1 if the requested element exists in the array, else it returns 0.
//Args: int* a-> Pointer to the array
//      begin ->    Helps us find the beginning position of the search algorithm
//      end->        Helps us find the ending position of the algorithm
//      elem->        The wanted element
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
