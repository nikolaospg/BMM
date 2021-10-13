#ifndef BMM_H
#define BMM_H

#include <sparse.h>
#include <block.h>

#include <omp.h>

#define AF_THRESH 3.5

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
    CSCMatrixfree(X);
    free(X);
}

/**
 * serial block BMM
 */
CSCMatrix* bmm_bp(CSCMatrix* A, CSCMatrix* B, int b) {
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

CSCMatrix* bmm_bpf(CSCMatrix* A, CSCMatrix* B, CSCMatrix* F, int b) {
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
    //#pragma omp parallel for schedule(static) collapse(2)
    for (int p = 0; p < n_b; p++) {
        for (int q = 0; q < n_b; q++) {
            int blockc_idx = p*n_b+q;
            // Cp,q = 0
            blockc[blockc_idx] = malloc(sizeof(CSCMatrix));
            blockc[blockc_idx]->n = b;
            blockc[blockc_idx]->col_ptr = calloc(b+1, sizeof(int));
            blockc[blockc_idx]->row_idx = NULL; // NULL can be passed to realloc later
            block_bmmf(blocka, blockb, blockf, &blockc[blockc_idx], p, q);

            //#pragma omp atomic
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
    } else if (strcmp(method, "bp") == 0) {
        if (filter) {
            return bmm_bpf(A, B, F, b);
        } else {
            return bmm_bp(A, B, b);
        }
    } else if (strcmp(method, "ss") == 0) {
        if (filter) {
            return bmm_ssf(A, B, F);
        } else {
            return bmm_ss(A, B);
        }
    } else if (strcmp(method, "sp") == 0) {
        if (filter) {
            printf("BMM SMMP parallel version not supported for filtered BMM. Use cp instead\n");
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