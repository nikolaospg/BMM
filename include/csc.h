#ifndef CSC_H
#define CSC_H

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#include "utils.h"
#include "mmio.h"

#include <omp.h>

#define AF_THRESH 3.5

/** 
 * Defining a Struct, which resembles the CSC data structure 
 * The matrices are boolean so we don't store the values in any array.
 **/
typedef struct{            
    int* row_idx;     // contains the row indices of each nonzero element
    int* col_ptr;     // contains the indices of row_idx in which there are nonzero elements for the requested column
    int n;            //number of columns (=number of rows because the matrix is square)
}CSCMatrix;

typedef struct{            
    int* col_idx;
    int* row_ptr;
    int n;
} CSRMatrix;

typedef struct{
    int nb;     //The number of blocks
    int b;      //The b of the blocks
    int n;      //The dimension of the original matrix
    int nnz;    //The number of non zero elements

    int* row_idx_combined;      //This array contains the row_idx of all the blocks, row_idx_combined=[row_idx1, row_idx2, row_idx3 ...] for every block.
    int* col_ptr_combined;      //This array contains the col_ptr of all the blocks, col_ptr_combined=[col_ptr1, col_ptr2, col_ptr3 ...] for every block.

    int* row_idx_indices;       //This helps us properly index the row_idx_combined array. If we have a block with block_row and block_col, 
    //then by using int* my_pointer= row_idx_combined + row_idx_indices[block_row*nb + block_col] we get a pointer to the appropriate place to access the row_idx of the block.

    int* col_ptr_indices;       //This helps us properly index the col_ptr_combined array. If we have a block with block_row and block_col, 
    //then by using int* my_pointer= col_ptr_combined + col_ptr_indices[block_row*nb + block_col] we get a pointer to the appropriate place to access the col_ptr of the block.


}CSCMatrixBlocked;

CSCMatrixBlocked* create_empty_block_matrix(int n, int b, int nnz) {
    CSCMatrixBlocked* ret = malloc(sizeof(CSCMatrixBlocked));
    int nb = ceil(n/((double) b));
    ret->nb = nb;
    ret->b = b;
    ret->n = n;
    ret->nnz = nnz;
    ret->row_idx_combined = malloc(nnz*sizeof(int));             
    ret->col_ptr_combined = malloc(nb*nb*(b+1)*sizeof(int));       
    ret->row_idx_indices = malloc((nb*nb+1)*sizeof(int));
    ret->col_ptr_indices = malloc((nb*nb+1)*sizeof(int));
    return ret;
}

/**
 * Frees the memory allocated for the row_idx and col_ptr of a specific CSCMatrix structure.
 **/
void CSCMatrixfree(CSCMatrix* arg){
    free(arg->col_ptr);
    free(arg->row_idx);
}

void CSRMatrixfree(CSRMatrix* arg){
    free(arg->row_ptr);
    free(arg->col_idx);
}

CSCMatrix* csc_copy(CSCMatrix* m) {
    CSCMatrix* ret = malloc(sizeof(CSCMatrix));
    ret->n = m->n;
    ret->col_ptr = malloc((m->n+1)*sizeof(int));
    ret->row_idx = malloc(m->col_ptr[m->n]*sizeof(int));
    memcpy(ret->col_ptr, m->col_ptr, (m->n+1)*sizeof(int));
    memcpy(ret->row_idx, m->row_idx, m->col_ptr[m->n]*sizeof(int));
    return ret;
}

/**
 * Compares two matrices for equivalency. It is not guaranteed that row indices are in ascending order.
 * Exits in case of failure with an error message.
 */
void csc_validate(CSCMatrix* A, CSCMatrix* B) {
    if (A->n != B->n) {
        printf("Error, matrices do not have equal length %d!=%d\n", A->n, B->n);
        exit(-1);
    }
    for (int i = 0; i <= A->n; i++) {
        if (A->col_ptr[i] != B->col_ptr[i]) {
            printf("Error, matrices do not have equal %d-th col_ptr elem %d!=%d\n", i, A->col_ptr[i], B->col_ptr[i]);
            exit(-1);
        }
    }

    int* mask = calloc(A->n, sizeof(int)); // contains whether a row_idx exists in the other's row_idx list
    for (int i = 0; i < A->n; i++) {
        // for each column, check if every element of B is contained in A and vice-versa
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) {
            int j = A->row_idx[jj];
            mask[j] = 1;
        }
        for (int jj = B->col_ptr[i]; jj < B->col_ptr[i+1]; jj++) {
            int j = B->row_idx[jj];
            if (!mask[j]) {
                printf("Error, second matrix does not contain row_idx %d in col %d\n", j, i);
                exit(-1);
            }
        }
        // clear mask and vice-versa
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) mask[A->row_idx[jj]] = 0;
        for (int jj = B->col_ptr[i]; jj < B->col_ptr[i+1]; jj++) {
            int j = B->row_idx[jj];
            mask[j] = 1;
        }
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) {
            int j = A->row_idx[jj];
            if (!mask[j]) {
                printf("Error, first matrix does not contain row_idx %d in col %d\n", j, i);
                exit(-1);
            }
        }
        // clear mask for next iteration
        for (int jj = B->col_ptr[i]; jj < B->col_ptr[i+1]; jj++) mask[B->row_idx[jj]] = 0;
    }
    free(mask);
}

/**
 * Converts CSC matrix to CSR
 * https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/csr.h#L418
 * O(nnz(A) + max(n_row,n_col))
 */
CSRMatrix* csc2csr(CSCMatrix* m) {
    int nnz = m->col_ptr[m->n];
    CSRMatrix* ret = malloc(sizeof(CSRMatrix));
    ret->n = m->n;
    ret->row_ptr = calloc(ret->n+1, sizeof(int));
    ret->col_idx = malloc(nnz * sizeof(int));

    for (int i = 0; i < nnz; i++)
        ret->row_ptr[m->row_idx[i]]++;

    //cumsum the nnz
    for(int row = 0, cumsum = 0; row < ret->n; row++){
        int temp  = ret->row_ptr[row];
        ret->row_ptr[row] = cumsum;
        cumsum += temp;
    }
    ret->row_ptr[ret->n] = nnz;

    for(int col = 0; col < m->n; col++){
        for(int jj = m->col_ptr[col]; jj < m->col_ptr[col+1]; jj++){
            int row  = m->row_idx[jj];
            int dest = ret->row_ptr[row];

            ret->col_idx[dest] = col;

            ret->row_ptr[row]++;
        }
    }

    for(int row = 0, last = 0; row <= ret->n; row++){
        int temp  = ret->row_ptr[row];
        ret->row_ptr[row] = last;
        last    = temp;
    }

    return ret;
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

//Used to free a CSCMatrixBlocked struct:
void CSCMatrixBlocked_free(CSCMatrixBlocked* A){
    free(A->row_idx_combined);
    free(A->col_ptr_combined);
    free(A->row_idx_indices);
    free(A->col_ptr_indices);
    free(A);
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

/*Function to convert a CSC structure to a row major array. Useful in validating the answers*/
int* CSR2array(CSRMatrix arg){
   /*Initialising useful variables and allocatin memory*/
   int n=arg.n; 
   int* ret=(int*)calloc(n*n,sizeof(int));
   if(ret==NULL){
       printf("could not allocate memory in CSC_convert function, exiting\n");
       exit(-1);
   }

   for(int i=0; i<n; i++){
       for(int jj=arg.row_ptr[i]; jj<arg.row_ptr[i+1]; jj++){
           int j=arg.col_idx[jj];
           ret[i*n+j]=1;
       }
   }
   return ret;
}

/* square CSC to CSR matrix conversion 
   https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L376
*/
CSRMatrix* CSC2CSR(CSCMatrix* m) {
    int nnz = m->col_ptr[m->n];
    CSRMatrix* ret = (CSRMatrix*) malloc(sizeof(CSRMatrix));
    ret->row_ptr = (int*) malloc((m->n+1)*sizeof(int));
    ret->col_idx = (int*) malloc(nnz*sizeof(int));
    ret->n = m->n;

    for (int i = 0; i < nnz; i++) {
        ret->row_ptr[m->row_idx[i]]++;
    }

    // nnz cumsum
    for (int i = 0, cumsum = 0; i < m->n; i++) {
        int tmp = ret->row_ptr[i];
        ret->row_ptr[i] = cumsum;
        cumsum += tmp;
    }
    ret->row_ptr[m->n] = nnz;

    // fill indices
    for(int col = 0; col < m->n; col++){
        for(int j = m->col_ptr[col]; j < m->col_ptr[col+1]; j++){
            int row  = m->row_idx[j];
            int dest = ret->row_ptr[row];

            ret->col_idx[dest] = col;
            ret->row_ptr[row]++;
        }
    }

    // restore row_ptr
    for(int row = 0, last = 0; row <= ret->n; row++){
        int tmp = ret->row_ptr[row];
        ret->row_ptr[row] = last;
        last = tmp;
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

void print_csr(CSRMatrix* m) {
    printf("N: %d\n", m->n);
    print_vector(m->row_ptr, m->n+1);
    print_vector(m->col_idx, m->row_ptr[m->n]);
    int* arr = CSR2array(*m);
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
    m->col_ptr=(int*)calloc((N+1), sizeof(int));
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
    return m;
}

//Fast version of the blocking. It returns a CSCMatrixBlocked struct
CSCMatrixBlocked* block_CSC(CSCMatrix* A, int b){
    //Useful variables
    int nb=ceil(A->n/((double) b));
    int n=A->n;
    int* col_ptr=A->col_ptr;
    int nnz=col_ptr[n];
    int* row_idx=A->row_idx;
    int current_block_row;
    int help_variable=(b+1)*nb;         //This is used on some products many times
    //Finished with useful variables


    //Now Allocating Memory For the arrays used in the CSCMatrixBlocked struct
    int* row_idx_combined=(int*)calloc(nnz, sizeof(int));             
    int* col_ptr_combined=(int*)calloc(nb*nb*(b+1), sizeof(int));       
    int* row_idx_indices= (int*)calloc(nb*nb +1, sizeof(int));
    int* col_ptr_indices= (int*)calloc(nb*nb +1, sizeof(int));

    if(row_idx_combined==NULL || col_ptr_combined ==NULL || row_idx_indices==NULL ||  col_ptr_indices==NULL){
        printf("Could not allocate memory in fast_block_CSC, exiting\n");
        exit(-1);
    }
    //Finished allocating memory


    //The following double for loop takes every non zero element of the CSCMatrix A, finds the block column and the block row of the element
    //and uses this info to initialise the values of the col_ptr_combined and the row_idx_indices arrays. These 2 now show the data in a non
    //cumulative way, in another loop later on they are changed so that they show the data in a cumulative way
    #pragma omp parallel for
    for (int current_block_col = 0; current_block_col < nb; current_block_col++) {
        for(int col=current_block_col*b; col<(current_block_col+1)*b; col++){
            if (col >= n) break;
            for (int index=col_ptr[col]; index<col_ptr[col+1]; index++ ){
                current_block_row=row_idx[index]/b;
                int access_index=current_block_row*help_variable + current_block_col +col;
                col_ptr_combined[access_index+1]++;
                row_idx_indices[current_block_row*nb + current_block_col +1]++;
            }
        }
    }

    //The following loop is used so that the col_ptr_combined shows the data in a cumulative way, as usual on the CSC scheme.
    #pragma omp parallel for schedule(static) collapse(2)
    for(int i=0; i<nb; i++){
        for(int j=0; j<nb; j++){
            int access_index=i*help_variable + j*(b+1) +1;
            for(int k=1; k<b+1; k++){
                col_ptr_combined[access_index]=col_ptr_combined[access_index] + col_ptr_combined[access_index-1];
                access_index++;
            }
        }
    }
    
    //The following for loop is used so that the row_idx_indices shows the data in a cumulative way, and also fills up the col_ptr_indices values
    int access_index=b+1;
    for(int i=1; i<nb*nb +1; i++){
        row_idx_indices[i]=row_idx_indices[i]+row_idx_indices[i-1];
        col_ptr_indices[i]=access_index;
        access_index=access_index+b+1;
    }

    //Accesing every non zero element to fill up the row_idx_combined array.
    int* count_array=(int*)calloc(nb*nb,sizeof(int));   //This is used so that I know how many elements I have stored on each block whenever needed.
    int offset;                                         //Another variable used to help me index the arrays
    for(int col=0; col<n; col++){
        int current_block_col=col/b;
        for (int index=col_ptr[col]; index<col_ptr[col+1]; index++ ){
            current_block_row=row_idx[index]/b;
            access_index=current_block_row*nb + current_block_col;
            offset=row_idx_indices[access_index] + count_array[access_index];
            row_idx_combined[offset]=row_idx[index]%b;
            count_array[access_index]++;
        }
    }
    free(count_array);
    //Finished with this too

    //Creating the CSCMatrixBlocked struct to be returned. For the explanation of what the members of the struct are, read the comments on the type definition
    CSCMatrixBlocked* ret=(CSCMatrixBlocked*)malloc(sizeof(CSCMatrixBlocked));
    ret->b=b;
    ret->nb=nb;
    ret->nnz=nnz;
    ret->n=n;
    ret->row_idx_combined=row_idx_combined;
    ret->row_idx_indices=row_idx_indices;
    ret->col_ptr_combined=col_ptr_combined;
    ret->col_ptr_indices=col_ptr_indices;
    //Finished creating the struct to be returned
    return ret;
}

//This works with a way equivalent to fast_block_CSC, but creates a CSCMatrixBlocked struct which has the row_idx and the col_ptr of each block in a column major
//way. This is done because the B matrix might be needed to work this way.
CSCMatrixBlocked* block_CSC_cols(CSCMatrix* A, int b){
    //Useful variables
    int nb=ceil(A->n/((double) b));
    int n=A->n;
    int* col_ptr=A->col_ptr;
    int nnz=col_ptr[n];
    int* row_idx=A->row_idx;
    int current_block_row;
    int current_row;
    int access_index;                   //I use this variable to help me access(index) some arrays with the proper way
    int help_variable=(b+1)*nb;         //This is used on some products many times
    //Finished with useful variables


    //Now Allocating Memory For the arrays used in the CSCMatrixBlocked struct
    int* row_idx_combined=(int*)calloc(nnz, sizeof(int));             
    int* col_ptr_combined=(int*)calloc(nb*nb*(b+1), sizeof(int));       
    int* row_idx_indices= (int*)calloc(nb*nb +1, sizeof(int));
    int* col_ptr_indices= (int*)calloc(nb*nb +1, sizeof(int));

    if(row_idx_combined==NULL || col_ptr_combined ==NULL || row_idx_indices==NULL ||  col_ptr_indices==NULL){
        printf("Could not allocate memory for block_CSC, exiting\n");
        exit(-1);
    }
    //Finished allocating memory


    //The following double for loop takes every non zero element of the CSCMatrix A, finds the block column and the block row of the element
    //and uses this info to initialise the values of the col_ptr_combined and the row_idx_indices arrays. These 2 now show the data in a non
    //cumulative way, in another loop later on they are changed so that they show the data in a cumulative way
    #pragma omp parallel for
    for (int current_block_col = 0; current_block_col < nb; current_block_col++) {
        for(int col=current_block_col*b; col<min((current_block_col+1)*b,n); col++){
            //int current_block_col=col/b;
            for (int index=col_ptr[col]; index<col_ptr[col+1]; index++ ){
                current_row=row_idx[index];
                current_block_row=current_row/b;
                access_index=current_block_col*help_variable + current_block_row*(b+1) +col%b;

                col_ptr_combined[access_index+1]++;
                row_idx_indices[current_block_col*nb + current_block_row +1]++;
            }
        }
    }
    //The following loop is used so that the col_ptr_combined shows the data in a cumulative way, as usual on the CSC scheme.
    #pragma omp parallel for schedule(static) collapse(2)
    for(int i=0; i<nb; i++){
        for(int j=0; j<nb; j++){
            access_index=i*help_variable + j*(b+1) +1;
            for(int k=1; k<b+1; k++){
                col_ptr_combined[access_index]=col_ptr_combined[access_index] + col_ptr_combined[access_index-1];
                access_index++;
            }
        }
    }
    
    //The following for loop is used so that the row_idx_indices shows the data in a cumulative way, and also fills up the col_ptr_indices values
    access_index=b+1;
    for(int i=1; i<nb*nb +1; i++){
        row_idx_indices[i]=row_idx_indices[i]+row_idx_indices[i-1];
        col_ptr_indices[i]=access_index;
        access_index=access_index+b+1;
    }

    //Accesing every non zero element to fill up the row_idx_combined array.
    int* count_array=(int*)calloc(nb*nb,sizeof(int));   //This is used so that I know how many elements I have stored on each block whenever needed.
    int offset;                                         //Another variable used to help me index the arrays
    for(int col=0; col<n; col++){
        int current_block_col=col/b;
        for (int index=col_ptr[col]; index<col_ptr[col+1]; index++ ){
            current_row=row_idx[index];
            current_block_row=current_row/b;
            //access_index=current_block_row*nb + current_block_col;
            access_index=current_block_col*nb + current_block_row;
            offset=row_idx_indices[access_index] + count_array[access_index];
            row_idx_combined[offset]=current_row%b;
            count_array[access_index]++;
        }
    }
    free(count_array);
    //Finished with this too

    //Creating the CSCMatrixBlocked struct to be returned. For the explanation of what the members of the struct are, read the comments on the type definition
    CSCMatrixBlocked* ret=(CSCMatrixBlocked*)malloc(sizeof(CSCMatrixBlocked));
    ret->b=b;
    ret->nb=nb;
    ret->nnz=nnz;
    ret->n=n;
    ret->row_idx_combined=row_idx_combined;
    ret->row_idx_indices=row_idx_indices;
    ret->col_ptr_combined=col_ptr_combined;
    ret->col_ptr_indices=col_ptr_indices;
    //Finished creating the struct to be returned
    return ret;
}

/**
 * Reconstructs a single CSC matrix from BCSC form
 */
CSCMatrix* unblock_CSC(CSCMatrixBlocked* blocked, int nb, int n) {
    int b=blocked->b;
    int nnz = blocked->nnz;
    int n_padded = nb*b;

    CSCMatrix* original = malloc(sizeof(CSCMatrix));
    original->n = n;
    original->col_ptr = calloc((n_padded+1), sizeof(int));
    original->col_ptr[0]=0;
    original->row_idx = malloc(nnz*sizeof(int));

    // first pass, original->col_ptr initially contains nnz per column
    #pragma omp parallel for schedule(static) collapse(2)
    for (int block_col = 0; block_col < nb; block_col++) { // for each block column
        for (int block_col_idx = 0; block_col_idx < b; block_col_idx++) { // for each column
            int col = block_col*b + block_col_idx;
            for (int block_row = 0; block_row < nb; block_row++) {
                int block_idx = block_row*nb + block_col;
                int block_col_ptr_off = blocked->col_ptr_indices[block_idx];
                int* block_col_ptr = blocked->col_ptr_combined + block_col_ptr_off;
                /*int block_row_idx_off = blocked->row_idx_indices[block_idx];
                int* block_row_idx = blocked->row_idx_combined + block_row_idx_off;*/
                
                original->col_ptr[col] += block_col_ptr[block_col_idx+1]-block_col_ptr[block_col_idx];
            }
        }
    }

    // cumsum to get final original->col_ptr
    for (int i = 0, cumsum = 0; i <= n_padded; i++) {
        int temp = cumsum;
        cumsum += original->col_ptr[i];
        original->col_ptr[i] = temp;
    }

    // second pass to compute original->row_idx
    #pragma omp parallel for schedule(static) collapse(2)
    for (int block_col = 0; block_col < nb; block_col++) { // for each block column
        for (int block_col_idx = 0; block_col_idx < b; block_col_idx++) { // for each column
            int col = block_col*b + block_col_idx;
            for (int block_row = 0; block_row < nb; block_row++) {
                int block_idx = block_row*nb + block_col;
                int block_col_ptr_off = blocked->col_ptr_indices[block_idx];
                int* block_col_ptr = blocked->col_ptr_combined + block_col_ptr_off;
                int block_row_idx_off = blocked->row_idx_indices[block_idx];
                int* block_row_idx = blocked->row_idx_combined + block_row_idx_off;
                
                for (int rrow = block_col_ptr[block_col_idx]; rrow < block_col_ptr[block_col_idx+1]; rrow++) {
                    int row = block_row*b + block_row_idx[rrow];
                    int dest = original->col_ptr[col];
                    original->row_idx[dest] = row;
                    original->col_ptr[col]++;
                }
            }
        }
    }
    // trim padding
    original->col_ptr = realloc(original->col_ptr, (n+1)*sizeof(int));

    // undo changes to original->col_ptr
    for(int col = 0, last = 0; col <= n; col++){
        int temp = original->col_ptr[col];
        original->col_ptr[col] = last;
        last = temp;
    }

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

/**
 * sparse vector multiplication. Assumes that indices are in ascending order
 */
int spvm(int* idx1, int begin1, int end1, int* idx2, int begin2, int end2) {
    int i = begin1;
    int j = begin2;
    
    while (i < end1 && j < end2) {
        if (idx1[i] == idx2[j])
            return 1;
        else if (idx1[i] < idx2[j])
            i++;
        else
            j++;
    }
    return 0;
}

/**
 * Computes CSC matrix element-wise or (C = A | B).
 * Assumes A->n=B->n=C->n and that C->col_ptr is correctly allocated with C->n elements
 * Assumes row_idx are in ascending order in each column. Do not perform in-line (with A=C or B=C)
 * 
 * 0(A.nnz+B.nnz)
 */
void spmor(CSCMatrix* A, CSCMatrix* B, CSCMatrix* C) {
    C->row_idx = realloc(C->row_idx, 
                (A->col_ptr[A->n]+B->col_ptr[B->n])*sizeof(int)); // allocate result upper bound for now
    int nnz = 0;
    C->col_ptr[0] = 0;

    for (int i = 0; i < C->n; i++) {
        int jja = A->col_ptr[i];
        int jjb = B->col_ptr[i];

        while(jja < A->col_ptr[i+1] && jjb < B->col_ptr[i+1]) {
            int ja = A->row_idx[jja];
            int jb = B->row_idx[jjb];

            if (ja < jb) {
                C->row_idx[nnz] = ja;
                jja++;
            } else if (jb < ja) {
                C->row_idx[nnz] = jb;
                jjb++;
            } else {
                C->row_idx[nnz] = ja;
                jja++;
                jjb++;
            }
            nnz++;
        }

        // Copy the remaining indices if left
        while (jja < A->col_ptr[i+1]) {
            int ja = A->row_idx[jja];
            C->row_idx[nnz] = ja;
            jja++;
            nnz++;
        }
        while (jjb < B->col_ptr[i+1]) {
            int jb = B->row_idx[jjb];
            C->row_idx[nnz] = jb;
            jjb++;
            nnz++;
        }
        C->col_ptr[i+1] = nnz;
    }

    C->row_idx = realloc(C->row_idx, nnz*sizeof(int));
}

/**
 * Computes CSC matrix element-wise or (C = A | B).
 * Assumes A->n=B->n=C->n and that C->col_ptr is correctly allocated with C->n elements
 * Does NOT assume row_idx are in ascending order in each column. Do not perform in-line (with A=C or B=C)
 * 
 * 0(A.nnz+B.nnz) using O(n) memory
 */
void spmor2(CSCMatrix* A, CSCMatrix* B, CSCMatrix* C) {
    C->row_idx = realloc(C->row_idx, 
                (A->col_ptr[A->n]+B->col_ptr[B->n])*sizeof(int)); // allocate result upper bound for now
    int nnz = 0;
    C->col_ptr[0] = 0;

    int* mask = calloc(C->n, sizeof(int)); // contains whether a row_idx exists in C's column or not

    for (int i = 0; i < C->n; i++) {
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) {
            int j = A->row_idx[jj];
            mask[j] = 1;
            C->row_idx[nnz] = j;
            nnz++;
        }
        for (int jj = B->col_ptr[i]; jj < B->col_ptr[i+1]; jj++) {
            int j = B->row_idx[jj];
            if (!mask[j]) {
                mask[j] = 1;
                C->row_idx[nnz] = j;
                nnz++;
            }
        }
        C->col_ptr[i+1] = nnz;
        // clear mask for next iteration
        for (int jj = C->col_ptr[i]; jj < C->col_ptr[i+1]; jj++) mask[C->row_idx[jj]] = 0;
    }
    free(mask);

    C->row_idx = realloc(C->row_idx, nnz*sizeof(int));
}
/*
 * Computes CSC matrix element-wise "and not" (C = !A & B).
 * Assumes A->n=B->n=C->n and that C->col_ptr is correctly allocated with C->n elements
 * Does NOT assume row_idx are in ascending order in each column. Do not perform in-line (with A=C or B=C)
 * 
 * 0(A.nnz+B.nnz) using O(n) memory
 */
void spmandnot(CSCMatrix* A, CSCMatrix* B, CSCMatrix* C) {
    C->row_idx = realloc(C->row_idx, 
                (A->col_ptr[A->n]+B->col_ptr[B->n])*sizeof(int)); // allocate result upper bound for now
    int nnz = 0;
    C->col_ptr[0] = 0;

    int* mask = calloc(C->n, sizeof(int)); // contains whether a row_idx exists in A's column or not

    for (int i = 0; i < C->n; i++) {
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) {
            int j = A->row_idx[jj];
            mask[j] = 1;
        }
        for (int jj = B->col_ptr[i]; jj < B->col_ptr[i+1]; jj++) {
            int j = B->row_idx[jj];
            if (!mask[j]) {
                C->row_idx[nnz] = j;
                nnz++;
            }
        }
        C->col_ptr[i+1] = nnz;
        // clear mask for next iteration
        for (int jj = A->col_ptr[i]; jj < A->col_ptr[i+1]; jj++) mask[A->row_idx[jj]] = 0;
    }
    free(mask);

    C->row_idx = realloc(C->row_idx, nnz*sizeof(int));
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

/**
 * BMM by converting A to CSR and doing CSR*CSC
 * Complexity: O(F.nnz*K) where K is the maximum nnz per column in A and B
 * Only viable with filtered BMM, otherwise O(F.n^2*K) is too expensive compared to bmm_ss
 * Uses OpenMP multithreading
 */
CSCMatrix* bmm_cpf(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F){
    CSRMatrix* a_csr = csc2csr(A);
    CSCMatrix* C = malloc(sizeof(CSCMatrix));
    C->col_ptr = calloc((B->n+1), sizeof(int));
    C->n = B->n;

    #pragma omp parallel for
    for (int j = 0; j < F->n; j++) {
        for (int ii = F->col_ptr[j]; ii < F->col_ptr[j+1]; ii++) {
            int i = F->row_idx[ii];
            int elem = spvm(a_csr->col_idx, a_csr->row_ptr[i], a_csr->row_ptr[i+1],
                            B->row_idx, B->col_ptr[j], B->col_ptr[j+1]);
            if (elem)
                C->col_ptr[j]++;
        }
    }

    for (int i = 0, cumsum = 0; i <= C->n; i++) {
        int temp = cumsum;
        cumsum += C->col_ptr[i];
        C->col_ptr[i] = temp;
    }
    C->row_idx = malloc(C->col_ptr[C->n]*sizeof(int));

    #pragma omp parallel for
    for (int j = 0; j < F->n; j++) {
        for (int ii = F->col_ptr[j]; ii < F->col_ptr[j+1]; ii++) {
            int i = F->row_idx[ii];
            int elem = spvm(a_csr->col_idx, a_csr->row_ptr[i], a_csr->row_ptr[i+1],
                             B->row_idx, B->col_ptr[j], B->col_ptr[j+1]);
            if (elem) {
                int dest = C->col_ptr[j];
                C->row_idx[dest] = i;
                C->col_ptr[j]++;
            }
        }
    }

    for(int col = 0, last = 0; col <= C->n; col++){
        int temp  = C->col_ptr[col];
        C->col_ptr[col] = last;
        last    = temp;
    }

    CSRMatrixfree(a_csr);
    free(a_csr);
    return C;
}

/**
 * BMM using SMMP algorithm for multiplication
 * https://www.researchgate.net/publication/2364309_Sparse_Matrix_Multiplication_Package_SMMP
 * Complexity: O(n*K^2) where K is the maximum nnz per column
 * Also uses O(n) space
 */
CSCMatrix* bmm_ss(CSCMatrix* A, CSCMatrix* B) {
    int n = B->n;
    CSCMatrix* C = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    C->n = n;
    C->col_ptr = (int*) malloc((n+1)*sizeof(int));

    // guess C's row_idx capacity
    int row_idx_capacity = B->n;
    C->row_idx = (int*) malloc(row_idx_capacity*sizeof(int));

    // used to not count multiple NZ per inner product index match
    int* mask = (int*) malloc(n*sizeof(int));
    for (int i = 0; i < n; i++) mask[i] = -1;
    C->col_ptr[0] = 0;

    int nnz = 0;
    for(int j = 0; j < n; j++){
        for(int kk = B->col_ptr[j]; kk < B->col_ptr[j+1]; kk++){
            int k = B->row_idx[kk];
            for(int ii = A->col_ptr[k]; ii < A->col_ptr[k+1]; ii++){
                int i = A->row_idx[ii];
                if(mask[i] != j){
                    mask[i] = j;
                    if (nnz >= row_idx_capacity) { // double row_idx capacity
                        row_idx_capacity *=2;
                        C->row_idx = realloc(C->row_idx, row_idx_capacity*sizeof(int));
                    }
                    C->row_idx[nnz] = i;
                    nnz++;
                }
            }
        }         
        C->col_ptr[j+1] = nnz;
    }
    free(mask);

    C->row_idx = realloc(C->row_idx, nnz*sizeof(int));
    return C;
}

/**
 * BMM using SMMP algorithm for multiplication, with filter
 * https://www.researchgate.net/publication/2364309_Sparse_Matrix_Multiplication_Package_SMMP
 * Complexity: O(n*K^3) where K is the maximum nnz per column
 * Also uses O(n) space
 */
CSCMatrix* bmm_ssf(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F) {
    int n = B->n;
    CSCMatrix* C = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    C->n = n;
    C->col_ptr = (int*) malloc((n+1)*sizeof(int));

    // guess C's row_idx capacity
    int row_idx_capacity = B->n;
    C->row_idx = (int*) malloc(row_idx_capacity*sizeof(int));

    // used to not count multiple NZ per inner product index match
    int* mask = (int*) malloc(n*sizeof(int));
    for (int i = 0; i < n; i++) mask[i] = -1;
    C->col_ptr[0] = 0;

    int nnz = 0;
    for(int j = 0; j < n; j++){
        for(int kk = B->col_ptr[j]; kk < B->col_ptr[j+1]; kk++){
            int k = B->row_idx[kk];
            for(int ii = A->col_ptr[k]; ii < A->col_ptr[k+1]; ii++){
                int i = A->row_idx[ii];
                for(int ll = F->col_ptr[j]; ll < F->col_ptr[j+1]; ll++){
                    int l = F->row_idx[ll];
                    if(l == i && mask[i] != j){
                        mask[i] = j;
                        if (nnz >= row_idx_capacity) { // double row_idx capacity
                            row_idx_capacity *=2;
                            C->row_idx = realloc(C->row_idx, row_idx_capacity*sizeof(int));
                        }
                        C->row_idx[nnz] = i;
                        nnz++;
                    }
                }
            }
        }         
        C->col_ptr[j+1] = nnz;
    }
    free(mask);

    C->row_idx = realloc(C->row_idx, nnz*sizeof(int));
    return C;
}

/**
 * BMM adaptive algorithm. Uses bmm_cpf for f.nnz/n>AF_THRESH and bmm_ssf otherwise
 */
CSCMatrix* bmm_af(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F) {
    int n = F->n;
    int nnz = F->col_ptr[n];
    double ratio = nnz / ((double) n);
    if (ratio > AF_THRESH)
        return bmm_cpf(A, B, F);
    else
        return bmm_ssf(A, B, F);
}

/**
 * BMM using SMMP algorithm, parallel version
 * https://www.researchgate.net/publication/2364309_Sparse_Matrix_Multiplication_Package_SMMP
 * Complexity: O(n*K^2) where K is the maximum nnz per column
 * Also uses O(n) space
 */
CSCMatrix* bmm_sp(CSCMatrix* A, CSCMatrix* B) {
    int n = B->n;
    CSCMatrix* C = (CSCMatrix*) malloc(sizeof(CSCMatrix));
    C->n = n;
    C->col_ptr = (int*) malloc((n+1)*sizeof(int));

    // used to not count multiple NZ per inner product index match
    int* mask = (int*) malloc(n*sizeof(int));
    #pragma omp parallel for
    for (int i = 0; i < n; i++) mask[i] = -1;
    C->col_ptr[0] = 0;

    int nnz = 0;
    for(int j = 0; j < n; j++){
        for(int kk = B->col_ptr[j]; kk < B->col_ptr[j+1]; kk++){
            int k = B->row_idx[kk];
            for(int ii = A->col_ptr[k]; ii < A->col_ptr[k+1]; ii++){
                int i = A->row_idx[ii];
                if(mask[i] != j){
                    mask[i] = j;
                    nnz++;
                }
            }
        }         
        C->col_ptr[j+1] = nnz;
    }
    C->row_idx = (int*) malloc(nnz*sizeof(int));

    // in the second pass, col_nnz holds the NNZ per column
    int* col_nnz = (int*) malloc(n*sizeof(int));
    // locks that protect mask array
    omp_lock_t* mask_lock = (omp_lock_t*) malloc(n*sizeof(omp_lock_t));
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        col_nnz[i] = C->col_ptr[i];
        mask[i] = -1;
        omp_init_lock(&mask_lock[i]);
    }
    //#pragma omp parallel for shared(mask, mask_lock) default(none)
    for (int j = 0; j < n; j++) {
        for(int kk = B->col_ptr[j]; kk < B->col_ptr[j+1]; kk++){
            int k = B->row_idx[kk];
            for(int ii = A->col_ptr[k]; ii < A->col_ptr[k+1]; ii++){
                int i = A->row_idx[ii];
                omp_set_lock(&mask_lock[i]);
                if (mask[i] != j) {
                    mask[i] = j;
                    C->row_idx[col_nnz[j]] = i;
                    col_nnz[j]++;
                }
                omp_unset_lock(&mask_lock[i]);
            }
        }         
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++) omp_destroy_lock(&mask_lock[i]);
    free(mask_lock);
    free(mask);
    free(col_nnz);
    return C;
}

/* Get submatrix of blocked in CSCMatrix* form */
CSCMatrix* get_matrix_block(CSCMatrixBlocked* blocked, int idx) {
    CSCMatrix* ret = malloc(sizeof(CSCMatrix));
    ret->col_ptr = blocked->col_ptr_combined + blocked->col_ptr_indices[idx];
    ret->row_idx = blocked->row_idx_combined + blocked->row_idx_indices[idx];
    ret->n = blocked->b;
    return ret;
}

/* Computes BMM of block Cpq from A_blocked and B_blocked*/
void block_bmm(CSCMatrixBlocked* A_blocked, CSCMatrixBlocked* B_blocked, CSCMatrix** Cpq, int p, int q) {
    int n_b = A_blocked->nb;
    for (int s = 0; s < n_b; s++) {
        int blocka_idx = p*n_b+s;
        int blockb_idx = s*n_b+q;
        CSCMatrix* ablock = get_matrix_block(A_blocked, blocka_idx);
        CSCMatrix* bblock = get_matrix_block(B_blocked, blockb_idx);

        CSCMatrix* AB = bmm_ss(ablock, bblock);
        CSCMatrix* tmp_blockc = csc_copy(*Cpq);
        spmor2(tmp_blockc, AB, *Cpq);
        CSCMatrixfree(AB);
        CSCMatrixfree(tmp_blockc);
        free(AB);
        free(tmp_blockc);

        free(ablock);
        free(bblock);
    }
}

/* Computes filtered BMM of block Cpq from A_blocked and B_blocked*/
void block_bmmf(CSCMatrixBlocked* A_blocked, CSCMatrixBlocked* B_blocked, CSCMatrixBlocked* F_blocked, CSCMatrix** Cpq, int p, int q) {
    int n_b = A_blocked->nb;
    int blockf_idx = p*n_b + q;
    CSCMatrix* fblock = get_matrix_block(F_blocked, blockf_idx);
    CSCMatrix* X = csc_copy(fblock);
    free(fblock);
    for (int s = 0; s < n_b; s++) {
        int blocka_idx = p*n_b+s;
        int blockb_idx = s*n_b+q;
        CSCMatrix* ablock = get_matrix_block(A_blocked, blocka_idx);
        CSCMatrix* bblock = get_matrix_block(B_blocked, blockb_idx);

        CSCMatrix* AB = bmm_af(ablock, bblock, X);
        CSCMatrix* tmp_blockc = csc_copy(*Cpq);
        spmor2(tmp_blockc, AB, *Cpq);
        CSCMatrixfree(AB);
        CSCMatrixfree(tmp_blockc);
        free(AB);
        free(tmp_blockc);

        CSCMatrix* tmp_x = csc_copy(X);
        spmandnot(*Cpq, tmp_x, X);
        CSCMatrixfree(tmp_x);
        free(tmp_x);

        free(ablock);
        free(bblock);
    }
}

/* converts a row major array of CSCMatrix* of size nb*nb to BSC format*/
CSCMatrixBlocked* blockcsc_tobsc(CSCMatrix** blocked, int nb, int n, int total_nnz) {
    CSCMatrixBlocked* bsc=(CSCMatrixBlocked*)malloc(sizeof(CSCMatrixBlocked));
    int b = blocked[0]->n;
    bsc->b=b;
    bsc->nb=nb;
    bsc->nnz=total_nnz;
    bsc->n=n;
    bsc->row_idx_combined=malloc(total_nnz*sizeof(int));
    bsc->row_idx_indices=malloc((nb*nb+1)*sizeof(int));
    bsc->col_ptr_combined=malloc(nb*nb*(b+1)*sizeof(int));
    bsc->col_ptr_indices=malloc((nb*nb+1)*sizeof(int));
    int nnz=0;
    for (int p = 0; p < nb; p++) {
        for (int q = 0; q < nb; q++) {
            int block_idx = p*nb+q;
            int block_nnz = blocked[block_idx]->col_ptr[b];
            memcpy(bsc->row_idx_combined+nnz, blocked[block_idx]->row_idx, block_nnz*sizeof(int));
            bsc->row_idx_indices[block_idx] = nnz;
            nnz += block_nnz;
            memcpy(bsc->col_ptr_combined+block_idx*(b+1), blocked[block_idx]->col_ptr, (b+1)*sizeof(int));
            bsc->col_ptr_indices[block_idx] = block_idx*(b+1);
        }
    }
    bsc->row_idx_indices[nb*nb] = nnz;
    bsc->col_ptr_indices[nb*nb] = nb*nb*(b+1);
    return bsc;
}

/**
 * serial block BMM
 */
CSCMatrix* bmm_bs(CSCMatrix* A, CSCMatrix* B, int b) {
    CSCMatrixBlocked* blocka = block_CSC(A, b);
    CSCMatrixBlocked* blockb = block_CSC(B, b);

    assert(A->n == B->n);

    int n_b=ceil(A->n/((double) b));
    CSCMatrix** blockc = (CSCMatrix**) malloc(n_b*n_b*sizeof(CSCMatrix*));   
    if(blockc == NULL){    
        printf("could not allocate memory in block_CSC function for the CSCMatrix** blocked, exiting\n");
        exit(-1);
    }
    int total_nnz = 0;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int p = 0; p < n_b; p++) {
        for (int q = 0; q < n_b; q++) {
            int blockc_idx = p*n_b+q;
            // Cp,q = 0
            blockc[blockc_idx] = malloc(sizeof(CSCMatrix));
            blockc[blockc_idx]->n = b;
            blockc[blockc_idx]->col_ptr = calloc(b+1, sizeof(int));
            blockc[blockc_idx]->row_idx = NULL; // NULL can be passed to realloc later
            block_bmm(blocka, blockb, &blockc[blockc_idx], p, q);
            
            #pragma omp atomic
            total_nnz += blockc[blockc_idx]->col_ptr[b];
        }
    }

    CSCMatrixBlocked* C_blocked = blockcsc_tobsc(blockc, n_b, A->n, total_nnz);
    for (int i = 0; i < n_b*n_b; i++) {
        CSCMatrixfree(blockc[i]);
        free(blockc[i]);
    }
    free(blockc);

    CSCMatrix* C = unblock_CSC(C_blocked, n_b, A->n);

    CSCMatrixBlocked_free(blocka);
    CSCMatrixBlocked_free(blockb);
    CSCMatrixBlocked_free(C_blocked);

    return C;
}

CSCMatrix* bmm_bsf(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F, int b) {
    CSCMatrixBlocked* blocka = block_CSC(A, b);
    CSCMatrixBlocked* blockb = block_CSC(B, b);
    CSCMatrixBlocked* blockf = block_CSC(F, b);

    assert(A->n == B->n);
    assert(A->n == F->n);

    int n_b=ceil(A->n/((double) b));
    CSCMatrix** blockc = (CSCMatrix**) malloc(n_b*n_b*sizeof(CSCMatrix*));   
    if(blockc == NULL){    
        printf("could not allocate memory in block_CSC function for the CSCMatrix** blocked, exiting\n");
        exit(-1);
    }
    int total_nnz = 0;
    #pragma omp parallel for schedule(static) collapse(2)
    for (int p = 0; p < n_b; p++) {
        for (int q = 0; q < n_b; q++) {
            int blockc_idx = p*n_b+q;
            // Cp,q = 0
            blockc[blockc_idx] = malloc(sizeof(CSCMatrix));
            blockc[blockc_idx]->n = b;
            blockc[blockc_idx]->col_ptr = calloc(b+1, sizeof(int));
            blockc[blockc_idx]->row_idx = NULL; // NULL can be passed to realloc later
            block_bmmf(blocka, blockb, blockf, &blockc[blockc_idx], p, q);

            #pragma omp atomic
            total_nnz += blockc[blockc_idx]->col_ptr[b];
        }
    }

    CSCMatrixBlocked* C_blocked = blockcsc_tobsc(blockc, n_b, A->n, total_nnz);
    for (int i = 0; i < n_b*n_b; i++) {
        CSCMatrixfree(blockc[i]);
        free(blockc[i]);
    }
    free(blockc);

    CSCMatrix* C = unblock_CSC(C_blocked, n_b, A->n);

    CSCMatrixBlocked_free(blocka);
    CSCMatrixBlocked_free(blockb);
    CSCMatrixBlocked_free(C_blocked);
    CSCMatrixBlocked_free(blockf);

    return C;
}

// interface for all BMM methods
CSCMatrix* bmm(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F, char* method, int filter, int b) {
    if (strcmp(method, "ds") == 0) {
        if (filter) {
            return bmm_dsf(A, B, F);
        } else {
            return bmm_ds(A, B);
        }
    } else if (strcmp(method, "bs") == 0) {
        if (filter) {
            return bmm_bsf(A, B, F, b);
        } else {
            return bmm_bs(A, B, b);
        }
    } else if (strcmp(method, "ss") == 0) {
        if (filter) {
            return bmm_ssf(A, B, F);
        } else {
            return bmm_ss(A, B);
        }
    } else if (strcmp(method, "sp") == 0) {
        if (filter) {
            printf("BMM SMMP parallel version not supported for non-filtered BMM. Use cp instead\n");
            exit(-1);
            return NULL;
        } else {
            return bmm_sp(A, B);
        }
    } else if (strcmp(method, "cp") == 0) {
        if (filter) {
            return bmm_cpf(A, B, F);
        } else {
            printf("BMM CSR*CSC not supported for non-filtered BMM. Use ss instead\n");
            exit(-1);
            return NULL;
        }
    } else if (strcmp(method, "a") == 0) {
        if (filter) {
            return bmm_af(A, B, F);
        } else {
            printf("BMM adaptive not supported for non-filtered BMM. Use ss instead\n");
            exit(-1);
            return NULL;
        }
    } else {
        fprintf(stderr, "Unknown method, %s, aborting\n", method);
        exit(-1);
    }
}

#endif
