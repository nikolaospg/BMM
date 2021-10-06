#include "csc.h"

//This script is used to demonstrate the use and correctness of the new bss functions. 


void partial_product(CSCMatrixBlocked* A, CSCMatrixBlocked* B, CSCMatrix* C, int s, int p, int q);

//Copies the content of block_c to temp_block_c
void csc_copy_new(CSCMatrix* block_c, CSCMatrix* temp_block_c) {

    temp_block_c->n = block_c->n;
    temp_block_c->row_idx = malloc(block_c->col_ptr[block_c->n]*sizeof(int));
    memcpy(temp_block_c->col_ptr, block_c->col_ptr, (block_c->n+1)*sizeof(int));
    memcpy(temp_block_c->row_idx, block_c->row_idx, block_c->col_ptr[block_c->n]*sizeof(int));
}


//This is the new ds function. Probably will be changed a little bit later on.
CSCMatrixBlocked* bmm_bs_new(CSCMatrix* A, CSCMatrix* B, int b) {

    //Blocking the A,B. This might change later (we might want to block them earlier)
    CSCMatrixBlocked* blocked_a=fast_block_CSC(A, b);
    CSCMatrixBlocked* blocked_b=fast_block_CSC_cols(B, b);
    assert(A->n == B->n);
    int n_b=ceil(A->n/((double) b));
    //Finished blocking

    //Initialising blocked_c and its arrays (the ones I already know)
    CSCMatrixBlocked* blocked_c=(CSCMatrixBlocked*)malloc(sizeof(CSCMatrixBlocked));
    if(blocked_c == NULL){    
        printf("could not allocate memory in bmm_bs function for the CSCMatrix* blocked_c, exiting\n");
        exit(-1);
    }       
    int* col_ptr_combined=(int*)calloc(n_b*n_b*(b+1), sizeof(int));       
    int* row_idx_indices= (int*)calloc(n_b*n_b +1, sizeof(int));
    int* col_ptr_indices= (int*)calloc(n_b*n_b +1, sizeof(int));
    if(col_ptr_combined ==NULL || row_idx_indices==NULL ||  col_ptr_indices==NULL){
        printf("Could not allocate memory for the arrays of blocked_c in bmm_bs, exiting\n");
        exit(-1);
    }
    int access_index=b+1;
    for(int i=1; i<n_b*n_b +1; i++){            //Filling up the values for col_ptr_indices
        col_ptr_indices[i]=access_index;
        access_index=access_index+b+1;
    }
    blocked_c->nb=n_b;
    blocked_c->b=b;
    blocked_c->n=A->n;      //The dimension of the unblocked matrices
    blocked_c->row_idx_indices=row_idx_indices;
    blocked_c->col_ptr_combined=col_ptr_combined;
    blocked_c->col_ptr_indices=col_ptr_indices;
    //Finished initialsing blocked_c

    
    //#pragma omp parallel for schedule(static) collapse(2)

    //initialising the AB struct, which holds the information for the A*B calculation. Probably will be put inside the for loop later on the parallel comp.
    CSCMatrix AB;
    AB.n=b;
    AB.col_ptr=(int*)calloc(b+1, sizeof(int));
    //Finished initialising the AB struct

    //initialising the temp_blockc struct, which holds the values of c. Probably will be put inside the for loop later on the parallel computation
    CSCMatrix temp_blockc;
    temp_blockc.n=b;
    temp_blockc.col_ptr=(int*)calloc(b+1, sizeof(int));
    //Finished initialising the temp_blockc struct

    //initialising the block_c struct, which implements the Cp,q we have in the equations
    CSCMatrix block_c;
    block_c.n=b;
    block_c.row_idx = NULL;
    //Finished initialising the temp_blockc struct

    //Initialising the two pointers that will help me merge the row_idxs of every Cp,q block
    int** row_idxs=(int**)malloc(n_b*n_b*sizeof(int*));
    //Finished initialising

    for (int p = 0; p < n_b; p++) {
        for (int q = 0; q < n_b; q++) {
            block_c.col_ptr=blocked_c->col_ptr_combined + blocked_c->col_ptr_indices[p*n_b +q];         //The col_ptr of block_c now points to the appropriate place on the col_ptr_combined
            block_c.row_idx = NULL;
            for (int s = 0; s < n_b; s++) {
                //In this case, both the A and B are the initial ones, on the MPI version they will be a subset of the original matrices and there will be a need of proper pointer arithmetic
                partial_product(blocked_a, blocked_b, &AB, s,p,q);      //With this command the AB values are changed
                csc_copy_new(&block_c, &temp_blockc);
                spmor2(&temp_blockc, &AB, &block_c);
                free(AB.row_idx);
                free(temp_blockc.row_idx);
            }
            row_idxs[p*n_b+q]=block_c.row_idx;
            blocked_c->row_idx_indices[p*n_b+q +1]=block_c.col_ptr[b];
        }
    }
    free(AB.col_ptr);
    free(temp_blockc.col_ptr);

    //Changing the row_idx_indices so that they show the data in a cumulative way
    for(int i=1; i<n_b*n_b +1; i++){
        blocked_c->row_idx_indices[i]=blocked_c->row_idx_indices[i]+blocked_c->row_idx_indices[i-1];
    }
    blocked_c->row_idx_indices[0]=0;

    //Merging the row_idxs, which contain the row_idxs for all the blocks, into blocked_c->row_idx_combined
    int nnz=row_idx_indices[n_b*n_b];
    blocked_c->nnz=nnz;
    blocked_c->row_idx_combined=(int*)malloc(nnz*sizeof(int));

    for(int i=0; i<n_b*n_b; i++){
        memcpy(blocked_c->row_idx_combined + row_idx_indices[i], row_idxs[i], (row_idx_indices[i+1]-row_idx_indices[i])*sizeof(int));
    }
    //Finished merging

    //Reconstructing
    CSCMatrixBlocked_free(blocked_a);               //If they are already blocked then we should not free them here
    CSCMatrixBlocked_free(blocked_b);
    for(int i=0; i<n_b*n_b; i++){
        free(row_idxs[i]);
    }
    free(row_idxs);

    
    return blocked_c;
}


//A is supposed to be a CSCMatrixBlocked with the pointers set in a way so that it represents a block row (or multiple ones), the memory for the 4 arrays is already allocated
//B is supposed to be a CSCMatrixBlocked with the pointers set in a way so that it represents a block col (or multiple ones), the memory for the 4 arrays is already allocated
//C is a block which in the end should represent AB. The col_ptr is already allocated. The only pointer which is not already allocated is the row_idx.
//p,q are the relative block_row and block_col. 
//The job of the function is to fill up the C->row_idx and C->col_ptr with the correct values
        //A0->col_ptr_combined=A2->col_ptr_combined+ (b+1)*nb*p;
        //A0->col_ptr_indices=A2->col_ptr_indices + nb*p;
        //A0->row_idx_indices=A2->row_idx_indices + nb*p;
        //A0->row_idx_combined=A2->row_idx_combined + A2->row_idx_indices[p*nb];
void partial_product(CSCMatrixBlocked* A, CSCMatrixBlocked* B, CSCMatrix* C, int s, int p, int q){

    //Useful Variables
    int nb=A->nb;
    int b=A->b;
    int* A_row_idx=A->row_idx_combined + A->row_idx_indices[s + p*nb] - A->row_idx_indices[0];
    int* A_col_ptr=A->col_ptr_combined + A->col_ptr_indices[s + p*nb] - A->col_ptr_indices[0];

    int* B_row_idx=B->row_idx_combined + B->row_idx_indices[s + q*nb] - B->row_idx_indices[0];
    int* B_col_ptr=B->col_ptr_combined + B->col_ptr_indices[s + q*nb] - B->col_ptr_indices[0];
    //Finished with the variables


    //guess C's row_idx capacity 
    int row_idx_capacity=b;
    C->row_idx = (int*) malloc(row_idx_capacity*sizeof(int));
    if(C->row_idx == NULL){
        printf("Could not allocate memory for C->row_idx on partial_product, exiting\n");
        exit(-1);
    }


    // used to not count multiple NZ per inner product index match
    int* mask = (int*) malloc(b*sizeof(int));
    if(mask == NULL){
        printf("Could not allocate memory for mask on partial_product, exiting\n");
        exit(-1);
    }

    for (int i = 0; i < b; i++) mask[i] = -1;
    C->col_ptr[0] = 0;

    int nnz = 0;
    for(int j = 0; j < b; j++){
        for(int kk = B_col_ptr[j]; kk < B_col_ptr[j+1]; kk++){
            int k = B_row_idx[kk];
            for(int ii = A_col_ptr[k]; ii < A_col_ptr[k+1]; ii++){
                int i = A_row_idx[ii];
                if(mask[i] != j){
                    mask[i] = j;
                    if (nnz >= row_idx_capacity) { // double row_idx capacity
                        row_idx_capacity *=2;
                        C->row_idx = realloc(C->row_idx, row_idx_capacity*sizeof(int));
                        if(C->row_idx == NULL){
                            printf("Could not reallocate memory for C->row_idx on partial_product, exiting\n");
                            exit(-1);
                        }
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
    if(C->row_idx == NULL){
        printf("Could not reallocate memory for C->row_idx on partial_product, exiting\n");
        exit(-1);
    }
}


int main(int argc, char* argv[]){
    srand(time(NULL));
    if(argc<4){
        printf("Give me the argument for the n, the sparsity, the method and the b value.\n");
        exit(-1);
    }
    int n=atoi(argv[1]);
    double sparsity;
    sscanf(argv[2], "%lf", &sparsity);
    int b_size=atoi(argv[3]);

    /*Initialising the useful matrices*/
    int* a=random_vector(n*n, sparsity);          //The two operand matrices, in 1D array form
    int* b=random_vector(n*n, sparsity);

    int* c=(int*)malloc(n*n*sizeof(int));       //The matrix from the classical multiplication

    //classical array multiplication
    int count=0;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            count=0;
            for(int k=0; k<n; k++){
                count=count+a[i*n+k]*b[j+k*n];
            }
            c[i*n+j]=count;
        }
    }
    //Doing the following normalisation, so that the result is boolean
    for(int i=0; i<n*n; i++){
        if(c[i]>0){
            c[i]=1;
        }
    }

    //Now doing the bss_new
    CSCMatrix A=array2CSC(a, n);
    CSCMatrix B=array2CSC(b,n);

    CSCMatrixBlocked* blocked_C=bmm_bs_new(&A, &B, b_size);
    CSCMatrix* C=fast_reconstruct_from_blocks(blocked_C);
    int* c_bmm=CSC2array(*C);            //Created the c_bmm matrix


    //Now comparing the two results
    for(int i=0; i<n*n; i++){
        //c[2]=-1;                  //Uncomment to deliberately create an error
        if(c_bmm[i]!= c[i]){
            printf("For i=%d c_bmm=%d c=%d on demo_bss, exiting\n",i,c_bmm[i], c[i]);
            exit(-1);
        }
    }


    free(a);
    free(b);
    free(c);

    CSCMatrixfree(&A);
    CSCMatrixfree(&B);
    CSCMatrixfree(C);
    CSCMatrixBlocked_free(blocked_C);
    free(c_bmm);
    free(C);

    return 0;
}
