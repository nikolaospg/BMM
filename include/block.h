#ifndef BLOCK_H
#define BLOCK_H

#include <sparse.h>

#include <omp.h>

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

//Used to free a CSCMatrixBlocked struct:
void CSCMatrixBlocked_free(CSCMatrixBlocked* A){
    free(A->row_idx_combined);
    free(A->col_ptr_combined);
    free(A->row_idx_indices);
    free(A->col_ptr_indices);
    free(A);
}

/* Get submatrix of blocked in CSCMatrix* form. */
CSCMatrix* get_matrix_block(CSCMatrixBlocked* blocked, int idx) {
    CSCMatrix* ret = malloc(sizeof(CSCMatrix));
    ret->col_ptr = blocked->col_ptr_combined + blocked->col_ptr_indices[idx];
    ret->row_idx = blocked->row_idx_combined + blocked->row_idx_indices[idx];
    ret->n = blocked->b;
    return ret;
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

/* converts a row major array of CSCMatrix* of size nb*nb to CSCMatricBlocked* */
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

#endif