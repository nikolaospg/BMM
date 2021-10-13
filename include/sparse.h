#ifndef SPARSE_H
#define SPARSE_H

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

typedef struct{            
    int* col_idx;
    int* row_ptr;
    int n;
} CSRMatrix;

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
 * compressed vector multiplication. Assumes that indices are in ascending order
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

#endif
