#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>



/** 
 * Defining a Struct, which resembles the CSC data structure 
 * The matrices are boolean so we don't store the values in any array.
 **/

typedef struct{            
    int* row_vector;     //array containing the row indices of each nonzero element
    int* col_vector;     //array containing the index of the elements which start a column of the sparse matrix
    int n;              //number of columns (=number of rows because the matrix is square)
}CSCArray;


/**
 * Function that frees the memory allocated for the row_vector and col_vector of a specific CSCArray structure.
 **/
void CSCArrayfree(CSCArray* arg){
    free(arg->col_vector);
    free(arg->row_vector);
}


/** 
 * Function to print a vector(row major fashion)
 * */
void print_vector(void* vector, int size){
    int* vector_temp=(int*)vector;
    for(int i=0; i<size; i++){
        printf("%d  ",vector_temp[i]);
    }
    printf("\n\n");

}



//Prints a vector as 2D matrix (row-major)
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




/**
 * Function that converts a row major array to CSC... Probably won't be included in the final project
 * inputs:init_array -> The initial array (row_major)
 *        n    ->  The number of rows/cols
 * outputs: The CSCArray structure
 * */
CSCArray array2CSC(int* init_array, int n){
    /*Allocating the memory for the structure members*/
    int* row_vector=(int*)malloc(n*sizeof(int));        //Initialising with this size to be sure, later on we will be reallocating. These are because I dont know the number of non zero elements a priory.
    if(row_vector==NULL){
        printf("could not allocate memory in array2CSC function, exiting\n");
        exit(-1);
    }
    int* col_vector=(int*)malloc((n+1)*sizeof(int));
    if(col_vector==NULL){
        printf("could not allocate memory in array2CSC function, exiting\n");
        exit(-1);
    }
    col_vector[0]=0;
    /*Finished with the memory allocation*/

    /*Double for loop to search every element of the initial array. Each i iteration corresponds to one column, each j iter. corresponds to one row.*/
    int row_vector_index=0; 
    int col_non_zero;       //The number of non zero elements of this current column
    for (int i=0; i<n; i++){
        col_non_zero=0;
        for(int j=0; j<n; j++){
           if(init_array[j*n+i]!=0){        //If the init_array was column major, then instead of j*n+i I would have had i*n+j
               row_vector[row_vector_index]=j;
               row_vector_index++;
               col_non_zero++;
           } 
        }
        col_vector[i+1]=col_vector[i]+col_non_zero;                             //Filling the col_vector with the correct values
        row_vector=(int*)realloc(row_vector, (col_vector[i+1]+n)*sizeof(int));      //Each time I finish with a line, I reallocate with the number of elements up until then plus n
    }
    /*Finished with the double for loop*/
    
    row_vector=(int*)realloc(row_vector, (col_vector[n])*sizeof(int));      //Final reallocation with the proper size
    CSCArray ret={row_vector, col_vector,n};
    return ret;
}



/*Function to convert a CSC structure to a row major array. Useful in validating the answers*/
int* CSC_convert(CSCArray arg){


   /*Initialising useful variables and allocatin memory*/
   int n=arg.n; 
   int* ret=(int*)calloc(n*n,sizeof(int));
   if(ret==NULL){
       printf("could not allocate memory in CSC_convert function, exiting\n");
       exit(-1);
   }
   int current_non_zeros;           //The number of non zero elements that a specific column will have.
   int row_count=0;                 //Counter of the rows in the row_vector of the structure.
   int current_row;
   /*Finished with the variables*/

   /*In this double for loop I get the values for the array. Each i iteration corresponds to one column and each j one to a row*/
   for(int i=0; i<n; i++){
       current_non_zeros=arg.col_vector[i+1]-arg.col_vector[i];
       for(int j=0; j<current_non_zeros; j++){
           current_row=arg.row_vector[row_count];
           row_count++;
           ret[i+current_row*n]=1;
       }
   }
   return ret;
}



//Returns the index of the wanted element. If it does not exist in the array, it returns -1.
//Args: int* a-> Pointer to the array
//      initial_index ->    The index of the first element of the array (the user just passes 0), it is only different for the recursive calls inside the function.
//      size->          The number of elements of the array
//      wanted->        The wanted element
//      The reason we pass the initial_index argument is to be able to know the actual index of the element we get. The algorithm works by cutting the array and passing the parts of it
//      as arguments, so in the recursive calls the information about the initial index is lost.
int binary_search(int* a,int initial_index, int size, int wanted){
    if(size>0){
        int mid_index=size/2;
        if(a[mid_index]==wanted){
            return 1;
        }

        else if(a[mid_index]>wanted){
            return binary_search(a, initial_index, mid_index, wanted);      //I call the binary search for the left part of the array. The a  pointer and the initial index remain the same
        }

        else{
            return binary_search(a+mid_index+1,initial_index+ mid_index+1, size-mid_index-1, wanted);
        }
    }
    
    return 0;
}



/*Searches for the element (row,col) in the A array. If it exists (is equal to 1), it returns 1 , else it returns 0
 *  inputs: The CSCArray struct A
 *          int row->   The desired row index
 *          int col->   The desired col index
 *  returns: An binary integer (0,1), according to whether the search was successful or not
 *  The function works be calling the binary_search with the proper arguments
 * */
int search_in_A(CSCArray A, int row, int col){
    int help_index=A.col_vector[col];                   //Index that helps me access the elements in the A array
    int current_non_zero=A.col_vector[col+1]-A.col_vector[col];         //The non zero elements of the column
    int ret;
    ret=binary_search(A.row_vector+help_index,0, current_non_zero, row);
    return ret;
}



/*Function to compute the inner product of a row of the A Array (in CSC), with a col of the B array (boolean product)
 *  Inputs: CSCArray A  ->The two matrices, in CSC form
 *          CSCArray B
 *          int row->   The row of the A matrix that is part of the operation
 *          int col->   The col of the B matrix that is part of the operation
 *  Outputs: A binary integer, according to whether there was at least one common element
 * */
int inner_product(CSCArray A, CSCArray B, int row, int col){
    

    /*First creating some useful variables*/
    int current_non_zeros=B.col_vector[col+1]-B.col_vector[col];    //Used to tell me how many non zero elements are in the column (of the B array)
    int help_index=B.col_vector[col];                               //Used so I can easily access the elements of the row_vector
    int current_row;        
    int ret_value=0;
    /*Finished with the variables*/

    /*This loop searches for every element of the column (so the search is done based on the column of B).*/
    for(int i=0; i<current_non_zeros; i++){
        current_row=B.row_vector[help_index+i];
        if(search_in_A(A,row,current_row)==1){      //If We get even one 1(successful search), then we stop to avoid the unnecessary iterations
            ret_value=1;
            break;
        }
    }
    return ret_value;
}




/*function to compute the product between two arrays in csc form
 *  inputs: cscarray a->    the a matrix
 *          cscarray b->    the b matrix
 *          cscarray f->    the f matrix(mask)
 *  the output is the boolean multiplication This is the version with the F filter
 * */
CSCArray bmmfiltered(CSCArray A, CSCArray B, CSCArray F){

    /*First creating useful variables and allocating memory*/
    int n=A.n;
    int* product_col_vector=(int*)malloc((n+1)*sizeof(int));    //Allocating memory for the column vector of the structure to be returned
    if(product_col_vector==NULL){
        printf("Could not allocate memory for the col_vector in the product function, exiting\n");
        exit(-1);
    }
    product_col_vector[0]=0;
    int* product_row_vector=(int*)malloc(n*sizeof(int));        //We allocate with this size just for the beginning. Later on, we will be reallocating to make it larger. The size is not known a priory.
    if(product_row_vector==NULL){
        printf("Could not allocate memory for the row_vector in the product function, exiting\n");
        exit(-1);
    }

    int current_non_zeros;              //Counter for the amount of non zero elements in each column (for the product structure).
    int F_col_non_zeros;                //Tells me the amount of non zero elements in the column of the mask (we use it to know how many j iterations there will be)
    int row_help_index;                 //Index used to easily access the elements 
    int row_count=0;                    //I use this counter to easily manipulate the elements of the product row_vector
    int inner_prod; 
    /*Finished with the memory allocation and the variables*/

    /*Each i iteration corresponds to one column, and then each j to a row*/
    for(int i=0; i<n; i++){     //For one column   
        F_col_non_zeros=F.col_vector[i+1]-F.col_vector[i];
        row_help_index=F.col_vector[i];
        current_non_zeros=0;

        for(int j=0; j<F_col_non_zeros; j++){
            inner_prod=inner_product(A,B,F.row_vector[row_help_index+j],i);         //We find the procuct with the row that the mask allows(only when the mask has non zero values do we calculate the product)
            if(inner_prod==1){
                product_row_vector[row_count]=F.row_vector[row_help_index+j];       //If the inner_prod is ==1 then I store the element in the row_vector, and increase my counters
                row_count++;
                current_non_zeros++;
            }
        }
        product_col_vector[i+1]=product_col_vector[i]+current_non_zeros;            //I update the col_vector of the product with the amount of non zero elements I had
        product_row_vector=(int*)realloc(product_row_vector, (product_col_vector[i+1]+n)*sizeof(int));      //I reallocate with the total amount of non zero elements up until now, plus n(to be sure).
    }
    /*Finished doing the calculations*/

    product_row_vector=(int*)realloc(product_row_vector, (product_col_vector[n])*sizeof(int));      //I reallocate with the final size.
    CSCArray C={product_row_vector, product_col_vector, n};

    return C;

}



/*function to compute the product between two arrays in csc form
 *  inputs: cscarray a->    the a matrix
 *          cscarray b->    the b matrix
 *  the output is the boolean multiplication. In this version, there is no F filter
 * */

CSCArray bmm(CSCArray A, CSCArray B){

    /*First creating useful variables and allocating memory*/
    int n=A.n;
    
    int* product_col_vector=(int*)malloc((n+1)*sizeof(int));
    if(product_col_vector==NULL){
        printf("Could not allocate memory for the col_vector in the product function, exiting\n");
        exit(-1);
    }
    product_col_vector[0]=0;

    int* product_row_vector=(int*)malloc(n*sizeof(int));        //We allocate this size just to begin, it will later be reallocated
    if(product_row_vector==NULL){
        printf("Could not allocate memory for the row_vector in the product function, exiting\n");
        exit(-1);
    }
    int current_non_zeros;      //Counter for the amount of non zero elements of the current column (to be stored in the returned array)
    int row_count=0;            //to easily manipulate the elements of the row_vector of the product struct.
    int inner_prod;
    /*Finished with the variables and the memory*/


    /*Each i iteration corresponds to one column, and then each j to a row*/
    for(int i=0; i<n; i++){     //For one column
        current_non_zeros=0;
        for(int j=0; j<n; j++){
            inner_prod=inner_product(A,B,j,i);
            if(inner_prod==1){
                product_row_vector[row_count]=j;    //If the inner_prod is ==1 then I store the element in the row_vector, and increase my counters
                row_count++;
                current_non_zeros++;

            }
        }
        product_col_vector[i+1]=product_col_vector[i]+current_non_zeros;     //I update the col_vector of the product with the amount of non zero elements I had
        product_row_vector=(int*)realloc(product_row_vector, (product_col_vector[i+1]+n)*sizeof(int));          //I reallocate with the total amount of non zero elements up until now, plus n(to be sure).
    }

    product_row_vector=(int*)realloc(product_row_vector, (product_col_vector[n])*sizeof(int));  //I reallocate with the final size
    CSCArray C={product_row_vector, product_col_vector, n};

    return C;

}



