#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

#include "utils.h"
#include "mmio.h"


/** 
 * Defining a Struct, which resembles the CSC data structure 
 * The matrices are boolean so we don't store the values in any array.
 **/
typedef struct{            
    int* row_idx;     // contains the row indices of each nonzero element
    int* col_ptr;     // contains the indices of row_idx in which there are nonzero elements for the requested column
    int n;            //number of columns (=number of rows because the matrix is square)
}CSCMatrix;


/**
 * Frees the memory allocated for the row_idx and col_ptr of a specific CSCMatrix structure.
 **/
void CSCMatrixfree(CSCMatrix* arg){
    free(arg->col_ptr);
    free(arg->row_idx);
}

/**
 * Frees the memory allocated by a CSCMatrix**. nb is the dimension of the CSCMatrix** (such that nb*nb is the total amount of elements of the CSCMatrix**)
 **/
void ArrayCSCMatrixfree(CSCMatrix** arg, int nb){
    for(int i=0; i< nb*nb; i++){
        CSCMatrixfree(arg[i]);
        free(arg[i]);
    }
    free(arg);
}

/**
 * Function that converts a row major array to CSC
 * inputs:init_array -> The initial array (row_major)
 *        n    ->  The number of rows/cols
 * outputs: The CSCMatrix structure
 * */
CSCMatrix array2CSC(int* init_array, int n){
    /*Allocating the memory for the structure members*/
    int* row_idx=(int*)malloc(n*sizeof(int));        //Initialising with this size to be sure, later on we will be reallocating. These are because I dont know the number of non zero elements a priory.
    if(row_idx==NULL){
        printf("could not allocate memory in array2CSC function, exiting\n");
        exit(-1);
    }
    int* col_ptr=(int*)malloc((n+1)*sizeof(int));
    if(col_ptr==NULL){
        printf("could not allocate memory in array2CSC function, exiting\n");
        exit(-1);
    }
    col_ptr[0]=0;
    /*Finished with the memory allocation*/

    /*Double for loop to search every element of the initial array. Each i iteration corresponds to one column, each j iter. corresponds to one row.*/
    int row_idx_index=0; 
    int col_non_zero;       //The number of non zero elements of this current column
    for (int i=0; i<n; i++){
        col_non_zero=0;
        for(int j=0; j<n; j++){
           if(init_array[j*n+i]!=0){        //If the init_array was column major, then instead of j*n+i I would have had i*n+j
               row_idx[row_idx_index]=j;
               row_idx_index++;
               col_non_zero++;
           } 
        }
        col_ptr[i+1]=col_ptr[i]+col_non_zero;                             //Filling the col_ptr with the correct values
        row_idx=(int*)realloc(row_idx, (col_ptr[i+1]+n)*sizeof(int));      //Each time I finish with a line, I reallocate with the number of elements up until then plus n
    }
    /*Finished with the double for loop*/
    
    row_idx=(int*)realloc(row_idx, (col_ptr[n])*sizeof(int));      //Final reallocation with the proper size
    CSCMatrix ret={row_idx, col_ptr,n};
    return ret;
}

/*Function to convert a CSC structure to a row major array. Useful in validating the answers*/
int* CSC2array(CSCMatrix arg){
   /*Initialising useful variables and allocatin memory*/
   int n=arg.n; 
   int* ret=(int*)calloc(n*n,sizeof(int));
   if(ret==NULL){
       printf("could not allocate memory in CSC_convert function, exiting\n");
       exit(-1);
   }
   int current_non_zeros;           //The number of non zero elements that a specific column will have.
   int row_count=0;                 //Counter of the rows in the row_idx of the structure.
   int current_row;
   /*Finished with the variables*/

   /*In this double for loop I get the values for the array. Each i iteration corresponds to one column and each j one to a row*/
   for(int i=0; i<n; i++){
       current_non_zeros=arg.col_ptr[i+1]-arg.col_ptr[i];
       for(int j=0; j<current_non_zeros; j++){
           current_row=arg.row_idx[row_count];
           row_count++;
           ret[i+current_row*n]=1;
       }
   }
   return ret;
}

void print_csc(CSCMatrix* m) {
    printf("N: %d\n", m->n);
    print_vector(m->col_ptr, m->n+1);
    print_vector(m->row_idx, m->col_ptr[m->n]);
    int* arr = CSC2array(*m);
    print_vector2D(arr, m->n, m->n);
    free(arr);
}

/**
 * Creates a CSCMatrix from a matrix market file
 * inputs:init_array -> The initial array (row_major)
 *        n    ->  The number of rows/cols
 * outputs: The CSCMatrix structure
 * */
CSCMatrix* CSCfromMM(char* filename){
    printf("Parsing %s\n", filename);
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    int M, N, nz;
    MM_typecode t;
    int res;
    res = mm_read_banner(fp, &t);
    if (res) {
        printf("Error reading banner: %d\n", res);
        exit(res);
    }
    if (!mm_is_sparse(t)) {
        printf("Error, matrix is not sparse\n");
        exit(2);
    }
    int pattern = 1;
    if (!mm_is_pattern(t)) {
        printf("Warning: matrix is not pattern\n");
        pattern = 0;
    }
    /*if (!mm_is_symmetric(t)) {
        printf("Error, matrix is not symmetric\n");
        exit(4);
    }*/
    res = mm_read_mtx_crd_size(fp, &M, &N, &nz);
    if (res) {
        printf("Error reading size information: %d\n", res);
        exit(res);
    }
    if (M != N) {
        printf("Error, matrix is not square\n");
        exit(5);
    }

    int wrn_nzoo = 1;
    CSCMatrix* m = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    m->n = N;
    m->col_ptr=(int*)malloc((N+1)*sizeof(int));
    if(!m->col_ptr){
        printf("Error, could not allocate col_ptr for MM file\n");
        exit(7);
    }
    m->col_ptr[0] = 0;
    int* col_index = (int*) malloc(nz*sizeof(int));
    int* row_index = (int*) malloc(nz*sizeof(int));
    printf("Matrix market file %s valid, parsing...\n", filename);

    // matrix market does not store upper triangular part on symmetrix matrices so code below has to be
    // unavoidably shit to support symmetric/skew-symmetric/hermitian matrices which happen
    // to be almost all the large matrices
    int n_extra = 0;
    int is_symmetric = mm_is_symmetric(t) || mm_is_hermitian(t) || mm_is_skew(t);
    if (is_symmetric) {
        row_index = (int*) realloc(row_index, (2*nz+N)*sizeof(int));
        col_index = (int*) realloc(col_index, (2*nz+N)*sizeof(int));
    }
    for (int i=0; i<nz; i++)
    {
        if (!pattern) {
            double val;
            fscanf(fp, "%d %d %lf\n", &row_index[i+n_extra], &col_index[i+n_extra], &val);
            if (val != 0 && val != 1) {
                if (wrn_nzoo) {
                    printf("Warning: non-zero-or-one element %lf, will be interpreted as 1, ", val);
                    printf("this warning is only shown once\n");
                    wrn_nzoo = 0;
                }
            }
        } else {
            fscanf(fp, "%d %d\n", &row_index[i+n_extra], &col_index[i+n_extra]);
        }

        row_index[i+n_extra]--;  /* adjust from 1-based to 0-based */
        col_index[i+n_extra]--;
        m->col_ptr[col_index[i+n_extra]+1]++; // col_ptr temporarily holds the number of nz in the previous column

        if (is_symmetric && row_index[i+n_extra] != col_index[i+n_extra]) { // non-diagonal entry
            n_extra++;

            row_index[i+n_extra] = col_index[i+n_extra-1];
            col_index[i+n_extra] = row_index[i+n_extra-1];
            
            m->col_ptr[col_index[i+n_extra]+1]++;
        }
    }
    if (is_symmetric) {
        row_index = (int*) realloc(row_index, (nz+n_extra)*sizeof(int));
        col_index = (int*) realloc(col_index, (nz+n_extra)*sizeof(int));
    }

    printf("Matrix market file %s parsed successfully\n", filename);
    if (is_symmetric)
        printf("N: %d, lower triangular NZ: %d, total nz: %d\n", N, nz, nz+n_extra);
    else
        printf("N: %d, nz: %d\n", N, nz);

    // cumsum
    for (int i = 0; i < N; i++) {
        m->col_ptr[i+1] += m->col_ptr[i];
    }

    // move row indices so row_index are in column major order
    // m.col_ptr holds the position in g.row_index that the next row_index of the i-th column should be placed in
    m->row_idx = (int*) malloc((nz+n_extra)*sizeof(int));
    for (int i = 0; i < nz+n_extra; i++) {
         int dst = m->col_ptr[col_index[i]];
         m->row_idx[dst] = row_index[i];
         m->col_ptr[col_index[i]]++;
    }

    // undo changes in col_index
    int prev = 0;
    for (int i = 0; i < N; i++) {
        int tmp = m->col_ptr[i];
        m->col_ptr[i] = prev;
        prev = tmp;
    }

    assert(m->col_ptr[N] == nz+n_extra);

    free(col_index);
    free(row_index);
    fclose(fp);
    printf("Sparse matrix from file %s created successfully\n", filename);
    return m;
}

/**
 *  Creates a block CSC matrix
 *  inputs:
 *      CSCMatrix* A  -> The initial CSCMatirx
 *      int b->       The size of the new blocks. The initial matrix is partitioned into nbxnb blocks, where nb=n/b must be an integer. Also nb must be >b.
 *  returns:
 *      CSCMatrix** blocked -> This value points to an 1D array, with (nbxnb) CSCMatrix* values.
 * */
CSCMatrix** block_CSC(CSCMatrix* A, int b){
    if(A->n%b!=0){     //Demanding that nb=n/b is an integer
        printf("Warning, in block_CSC, the mod(n,b)=%d, which is !=0.\n",A->n%b);
    }
    int n_b=ceil(A->n/((double) b));
    printf("Creating block matrix with %dx%d blocks of size %dx%d\n", n_b, n_b, b, b);
    if(n_b<=b){      //Demanding that nb>b
        printf("Warning, the nb=%d, b=%d. nb is expected to be >b\n",n_b,b);
    }

    /*Now allocating memory for the CSCMatrix** (to be returned)*/
    CSCMatrix** blocked=(CSCMatrix**)malloc(n_b*n_b*sizeof(CSCMatrix*));   
    if(blocked==NULL){          
        printf("could not allocate memory in block_CSC function for the CSCMatrix** blocked, exiting\n");
        exit(-1);
    }    

    //In each iteration of the following for loop, We allocate memory for one of the nbxnb blocks. 
    int init_row_idx_capacity = A->col_ptr[A->n]/(n_b*n_b)+1; // expected capacity (+1 for non zero)
    int* row_idx_capacities = (int*) malloc(n_b*n_b*sizeof(int)); // contains the allocation size of row_idx for each block
    if(row_idx_capacities==NULL){
        printf("could not allocate memory in block_CSC function for row_idx capacities, exiting\n");
        exit(-1);
    }
    for (int i = 0; i < n_b*n_b; i++) row_idx_capacities[i] = init_row_idx_capacity;
    for (int i=0; i<n_b*n_b; i++){
        blocked[i]=(CSCMatrix*)malloc(sizeof(CSCMatrix));
        if(blocked[i]==NULL){
            printf("could not allocate memory in block_CSC function for the %dth element of the blocked array, exiting\n",i);
            exit(-1);
        }
        blocked[i]->n=b;          

        blocked[i]->col_ptr=(int*)malloc((b+1)*sizeof(int));    //I know that every block has b columns
        if(blocked[i]->col_ptr==NULL){
            printf("could not allocate memory in block_CSC function for the col_ptr of the %dth element of the blocked array, exiting\n",i);
            exit(-1);
        }
        blocked[i]->col_ptr[0]=0;
        blocked[i]->row_idx=(int*)malloc((row_idx_capacities[i])*sizeof(int));
        if(blocked[i]->row_idx==NULL){
            printf("could not allocate memory in block_CSC function for the row_idx of the %dth element of the blocked array, exiting\n",i);
            exit(-1);
        }
    }

    // fill blocks
    int n_padded = n_b*b;
    for (int i = 0; i < n_padded; i++) {
        int block_idx_x = i/b;
        int block_col_idx = i%b;
        for (int block_idx_y = 0; block_idx_y < n_b; block_idx_y++) {
            int block_idx = block_idx_y*n_b + block_idx_x; // block index in blocked (blocked is row-major)
            blocked[block_idx]->col_ptr[block_col_idx+1] = blocked[block_idx]->col_ptr[block_col_idx];
        }

        if (i < A->n) { // zero padding after A->n-th column
            for (int j = A->col_ptr[i]; j < A->col_ptr[i+1]; j++) {
                int row_idx = A->row_idx[j];
                int block_row_idx = row_idx%b;
                int block_idx_y = row_idx/b;
                int block_idx = block_idx_y*n_b + block_idx_x; // block index in blocked (blocked is row-major)

                int block_nnz = blocked[block_idx]->col_ptr[block_col_idx+1];
                int block_row_idx_capacity = row_idx_capacities[block_idx];
                if (block_nnz == block_row_idx_capacity) {
                    // expand row_idx capacity
                    row_idx_capacities[block_idx] *= 2;
                    blocked[block_idx]->row_idx = (int*) realloc(blocked[block_idx]->row_idx, row_idx_capacities[block_idx]*sizeof(int));
                }
                blocked[block_idx]->row_idx[block_nnz] = block_row_idx;
                blocked[block_idx]->col_ptr[block_col_idx+1]++;
            }
        }

        // after last column, reallocate blocks to just neccessary size
        if (block_col_idx == b-1) {
            for (int block_idx_y = 0; block_idx_y < n_b; block_idx_y++) {
                int block_idx = block_idx_y*n_b + block_idx_x; // block index in blocked (blocked is row-major)
                blocked[block_idx]->row_idx = (int*) realloc(blocked[block_idx]->row_idx, blocked[block_idx]->col_ptr[b]*sizeof(int));
                // row_idx_capacities[block_idx] = blocked[block_idx]->col_ptr[b]; // not neccessary, won't be used again
            }
        }
    }
    free(row_idx_capacities);
    
    return blocked;
}


/**
 *  The Following function helps us in the task of reconstructing one Matrix (CSC Form) from a blocked version of the same matrix
 *  It is used to verify that the process of blocking is done correctly
 *  inputs:
 *      1)CSCMatrix** blocked   -> the blocked version of the matrix
 *      2)int nb->      The nb parameter
 *  outputs:
 *      CSCMatrix* reconstructed-> The reconstructed CSC Matrix
 * */
CSCMatrix* reconstruct_from_blocks(CSCMatrix** blocked, int nb){

    /*Initialising useful variables*/
    int b=blocked[0]->n;
    int n=b*nb;
    int total_nz=0;

    //Calculating the total amount of non zero elements:
    for (int i=0; i<nb*nb; i++){
        total_nz=total_nz+blocked[i]->col_ptr[b];
    }

    //Initialising memory for the original matrix//
    CSCMatrix* original=(CSCMatrix*)malloc(sizeof(CSCMatrix));
    if(original==NULL){
        printf("On reconstruct_from_blocks could not allocate memory for original!, error\n");
        exit(-1);
    }

    original->col_ptr=(int*)malloc((n+1)*sizeof(int));
    if(original->col_ptr==NULL){
        printf("On reconstruct_from_blocks could not allocate memory for original->col_ptr!, error\n");
        exit(-1);
    }
    original->col_ptr[0]=0;


    original->row_idx=(int*)malloc((total_nz)*sizeof(int));
    if(original->row_idx==NULL){
        printf("On reconstruct_from_blocks could not allocate memory for original->col_ptr!, error\n");
        exit(-1);
    }

    original->n=n;
    //Finished initialising memory for the original matrix//


    CSCMatrix* current_block;
    int block_col;              //The BLOCK COLUMN where some specific block belongs to.
    int current_nz;         //The nz elements of a specific column on a specific block.
    int relative_col;       //A column of the non blocked matrix corresponds to a column of the blocks. We name his relative_col. It is given by relative_col=column%b.
    int current_relative_row;        //The element we just extracted (is a row on some block)
    int current_row;        //current_row= block_row*b + current_relative_row
    int row_count;            //counting how many elements I have found up until a specific point for one column. Again used for indexing.
    int column_offset;      //Helps us index the row_idx vector for the reconstructed matrix
    /*Finished initialising useful variables*/


    /*Now extracting the elements from the blocked Matrix and storing them in the reconstructed one */

    //Each col_index iteration of the following loop corresponds to one of the n columns of the reconstructed matrix
    for(int col_index=0; col_index<n; col_index++){
        column_offset=original->col_ptr[col_index];
        row_count=0;                            //This is set=0 for each new column we study

        //Each block row iteration of the following loop corresponds to one of the nb blocks in this specific column
        for(int block_row=0; block_row<nb; block_row++){
            block_col=col_index/b;                              //In which one of the nb BLOCK COLUMNS the current block belongs to.
            current_block=blocked[block_col+ block_row*nb];     //Using the current_block variable to point to the current block
            relative_col=col_index%b;
            current_nz=current_block->col_ptr[relative_col+1]- current_block->col_ptr[relative_col];        //The non zero elements in this column of the block (to be stored in the reconstructed)

            original->col_ptr[col_index+1]=original->col_ptr[col_index+1] +current_nz;      //Increasing the value of the col_index vector

            //Each i iteration of the following loop corresponds to one non zero element of the current column of the current block
            for(int i=0; i<current_nz; i++){
                current_relative_row=current_block->row_idx[current_block->col_ptr[relative_col]+i];
                current_row=block_row*b + current_relative_row;
                original->row_idx[column_offset+row_count]=current_row;
                row_count=row_count+1;


            }
        }
        original->col_ptr[col_index+1]=original->col_ptr[col_index]+row_count;

    }
    /*Finished with the reconstructing task*/
    
    return original;

}




/*Searches for the element (row,col) in the A array. If it exists (is equal to 1), it returns 1 , else it returns 0
 *  inputs: The CSCMatrix struct A
 *          int row->   The desired row index
 *          int col->   The desired col index
 *  returns: An binary integer (0,1), according to whether the search was successful or not
 *  The function works be calling the binary_search with the proper arguments
 * */
int CSC_read_elem(CSCMatrix* A, int row, int col){
    int col_start=A->col_ptr[col];
    int col_end=A->col_ptr[col+1];

    return binary_search(A->row_idx, col_start, col_end, row);
}

/*Function to compute the inner product of a row of the A Matrix (in CSC), with a col of the B Matrix (boolean product)
 *  Inputs: CSCMatrix A  ->The two matrices, in CSC form
 *          CSCMatrix B
 *          int row->   The row of the A matrix that is part of the operation
 *          int col->   The col of the B matrix that is part of the operation
 *  Outputs: A binary integer, according to whether there was at least one common element
 * */
int inner_product(CSCMatrix* A, CSCMatrix* B, int row, int col){
    /*First creating some useful variables*/
    int current_non_zeros=B->col_ptr[col+1]-B->col_ptr[col];    //Used to tell me how many non zero elements are in the column (of the B array)
    int help_index=B->col_ptr[col];                               //Used so I can easily access the elements of the row_idx
    int current_row;        
    int ret_value=0;
    /*Finished with the variables*/

    /*This loop searches for every element of the column (so the search is done based on the column of B).*/
    for(int i=0; i<current_non_zeros; i++){
        current_row=B->row_idx[help_index+i];
        if(CSC_read_elem(A,row,current_row)==1){      //If We get even one 1(successful search), then we stop to avoid the unnecessary iterations
            ret_value=1;
            break;
        }
    }
    return ret_value;
}


/*  BMM using definition, serial implementation, filtered
 *  inputs: CSCMatrix* A->    the a matrix
 *          CSCMatrix* B->    the b matrix
 *          CSCMatrix* F->    the f matrix(mask)
 *  the output is the boolean multiplication (CSCMatrix*)
 * */
CSCMatrix* bmm_dsf(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F){
    /*First creating useful variables and allocating memory*/
    int n=A->n;
    int* product_col_ptr=(int*)malloc((n+1)*sizeof(int));    //Allocating memory for the column vector of the structure to be returned
    if(product_col_ptr==NULL){
        printf("Could not allocate memory for the col_ptr in the product function, exiting\n");
        exit(-1);
    }
    product_col_ptr[0]=0;
    int* product_row_idx=(int*)malloc(n*sizeof(int));        //We allocate with this size just for the beginning. Later on, we will be reallocating to make it larger. The size is not known a priory.
    if(product_row_idx==NULL){
        printf("Could not allocate memory for the row_idx in the product function, exiting\n");
        exit(-1);
    }

    int current_non_zeros;              //Counter for the amount of non zero elements in each column (for the product structure).
    int F_col_non_zeros;                //Tells me the amount of non zero elements in the column of the mask (we use it to know how many j iterations there will be)
    int row_help_index;                 //Index used to easily access the elements 
    int row_count=0;                    //I use this counter to easily manipulate the elements of the product row_idx
    int inner_prod; 
    /*Finished with the memory allocation and the variables*/

    /*Each i iteration corresponds to one column, and then each j to a row*/
    for(int i=0; i<n; i++){     //For one column   
        F_col_non_zeros=F->col_ptr[i+1]-F->col_ptr[i];
        row_help_index=F->col_ptr[i];
        current_non_zeros=0;

        for(int j=0; j<F_col_non_zeros; j++){
            inner_prod=inner_product(A,B,F->row_idx[row_help_index+j],i);         //We find the procuct with the row that the mask allows(only when the mask has non zero values do we calculate the product)
            if(inner_prod==1){
                product_row_idx[row_count]=F->row_idx[row_help_index+j];       //If the inner_prod is ==1 then I store the element in the row_idx, and increase my counters
                row_count++;
                current_non_zeros++;
            }
        }
        product_col_ptr[i+1]=product_col_ptr[i]+current_non_zeros;            //I update the col_ptr of the product with the amount of non zero elements I had
        product_row_idx=(int*)realloc(product_row_idx, (product_col_ptr[i+1]+n)*sizeof(int));      //I reallocate with the total amount of non zero elements up until now, plus n(to be sure).
    }
    /*Finished doing the calculations*/

    product_row_idx=(int*)realloc(product_row_idx, (product_col_ptr[n])*sizeof(int));      //I reallocate with the final size.
    CSCMatrix* C = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    C->row_idx = product_row_idx;
    C->col_ptr = product_col_ptr;
    C->n = n;

    return C;
}

/*  BMM using definition, serial implementation
 *  inputs: CSCMatrix* A->    the a matrix
 *          CSCMatrix* B->    the b matrix
 *  the output is the boolean multiplication (CSCMatrix*). In this version, there is no F filter
 * */
CSCMatrix* bmm_ds(CSCMatrix* A, CSCMatrix* B){
    /*First creating useful variables and allocating memory*/
    int n=A->n;
    
    int* product_col_ptr=(int*)malloc((n+1)*sizeof(int));
    if(product_col_ptr==NULL){
        printf("Could not allocate memory for the col_ptr in the product function, exiting\n");
        exit(-1);
    }
    product_col_ptr[0]=0;

    int* product_row_idx=(int*)malloc(n*sizeof(int));        //We allocate this size just to begin, it will later be reallocated
    if(product_row_idx==NULL){
        printf("Could not allocate memory for the row_idx in the product function, exiting\n");
        exit(-1);
    }
    int current_non_zeros;      //Counter for the amount of non zero elements of the current column (to be stored in the returned array)
    int row_count=0;            //to easily manipulate the elements of the row_idx of the product struct.
    int inner_prod;
    /*Finished with the variables and the memory*/


    /*Each i iteration corresponds to one column, and then each j to a row*/
    for(int i=0; i<n; i++){     //For one column
        current_non_zeros=0;
        for(int j=0; j<n; j++){
            inner_prod=inner_product(A,B,j,i);
            if(inner_prod==1){
                product_row_idx[row_count]=j;    //If the inner_prod is ==1 then I store the element in the row_idx, and increase my counters
                row_count++;
                current_non_zeros++;

            }
        }
        /*double percent = 100.0*((double)i)/n;
        printf("done a column %lf\n", percent);*/
        product_col_ptr[i+1]=product_col_ptr[i]+current_non_zeros;     //I update the col_ptr of the product with the amount of non zero elements I had
        product_row_idx=(int*)realloc(product_row_idx, (product_col_ptr[i+1]+n)*sizeof(int));          //I reallocate with the total amount of non zero elements up until now, plus n(to be sure).
    }

    product_row_idx=(int*)realloc(product_row_idx, (product_col_ptr[n])*sizeof(int));  //I reallocate with the final size
    CSCMatrix* C = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    C->row_idx = product_row_idx;
    C->col_ptr = product_col_ptr;
    C->n = n;

    return C;
}

// interface for all BMM methods
CSCMatrix* bmm(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F, char* method, int filter) {
    if (strcmp(method, "ds") == 0) {
        if (filter) {
            return bmm_dsf(A, B, F);
        } else {
            return bmm_ds(A, B);
        }
    } else {
        fprintf(stderr, "Unknown method, %s, aborting\n", method);
        exit(-1);
    }
}
