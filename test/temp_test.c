#include "csc.h"

//TEMPORARY SCRIPT. TO EVALA GIA NA ELENKSEIS KAI ESI AN DOULEVEI KALA O NEOS ALGORITHMOS GIA BLOCKING KAI GIA NA PSILODEIS LIGO POS LEITOURGEI
//DEN EXEI KAI IDIAITERA MEGALES DIAFORES SE SXESI ME PRIN!!!!

//Function which takes the 2 different blocked matrices and compares every single block, so that we know if the new implementation is ok
void compare_block_functions(CSCMatrixBlocked* A1, CSCMatrix** A2){
    
    //First printing the scalar values of A1, so that we can check whether they are correct or not
    printf("nb=%d b=%d n=%d nnz=%d \n",A1->nb, A1->b, A1->n, A1->nnz);

    //Initialising useful variables
    int nb=A1->nb;
    int b=A1->b;
    int current_non_zeros;      //The non zero elements of some block we will be accessing.
    CSCMatrix* temp2;       //Helps us get the proper A2 block each time
    int* temp1_row_idx;     //This pointer helps me get the proper row_idx subarray. It does it by doing THE PROPER pointer arithmetic, which is explained on the definition of the CSCMatrixBlocked struct 
    int* temp1_col_ptr;     //This pointer helps me get the proper col_ptr subarray, again with the proper pointer arithmetic.
    //Finished initialising useful variables
    

    //Each i iteration of the following for loop corresponds to one block row
    for(int i=0; i<nb; i++){
        //Each j iteration of the following for loop corresponds to one block column
        //This way, Each (i,j) corresponds to one block 
        for(int j=0; j<nb; j++){
            temp2=A2[i*nb +j];              //Getting the A2 block corresponding to this (i,j)

            temp1_row_idx=A1->row_idx_combined + A1->row_idx_indices[i*nb + j];             //USING POINTER ARITHMETIC TO POINT TO THE PROPER PLACE - THE PROPER row_idx SUBARRAY
            temp1_col_ptr=A1->col_ptr_combined + A1->col_ptr_indices[i*nb + j];             //USING POINTER ARITHMETIC TO POINT TO THE PROPER PLACE - THE PROPER col_ptr SUBARRAY

            //Each iteration of the following loop compares the col_ptr of this pair of blocks
            for(int col_index=0; col_index<b+1; col_index++){
                //temp1_col_ptr[col_index]=-1;                  //Uncomment to deliberately create an error
                if( temp1_col_ptr[col_index] != temp2->col_ptr[col_index]){
                    printf("ERROR i=%d j=%d col_index=%d, temp1=%d temp2=%d\n",i,j,col_index, temp1_col_ptr[col_index], temp2->col_ptr[col_index]);
                    exit(-1);
                }
            }

            
            current_non_zeros=temp2->col_ptr[b];
            for(int row_index=0; row_index<current_non_zeros; row_index++){
                //temp1_row_idx[row_index]=-1;                  //Uncomment to deliberately create an error
                if( temp1_row_idx[row_index] != temp2->row_idx[row_index]){
                    printf("ERROR i=%d j=%d row_index=%d, temp1=%d temp2=%d\n",i,j,row_index, temp1_col_ptr[row_index], temp2->col_ptr[row_index]);
                    exit(-1);
                }
            }
        }
    }

}

//This function compares the two reconstruct functions.
void compare_reconstruct_functions(CSCMatrix* fast, CSCMatrix* first){

    //Comparing the n
    if(fast->n!=first->n){
        printf("Different n! exiting\n");
        exit(-1);
    }

    
    for(int i=0; i<first->col_ptr[fast->n]; i++){
        //fast->row_idx[i]=-1;        //Uncomment to deliberately create an error
        if(fast->row_idx[i] != first->row_idx[i]){
            printf("Different row_idx[%d] on compare_reconstruct_functions! exiting\n",i);
            exit(-1);
        }
    }

    for(int i=0; i<first->n+1; i++){
        //fast->col_ptr[i]=-1;        //Uncomment to deliberately create an error
        if(fast->col_ptr[i] != first->col_ptr[i]){
            printf("Different col_ptr[%d] on compare_reconstruct_functions! exiting\n",i);
            exit(-1);
        }
    }


}

//This function compares the fast_block_CSC with the fast_block_CSC_cols function.
void compare_cols_rows(CSCMatrixBlocked* row_major, CSCMatrixBlocked* col_major){

    //Initialising useful variables
    int b=row_major->b;
    int current_non_zeros;
    int* row_major_row_idx;
    int* row_major_col_ptr;

    int* col_major_row_idx;
    int* col_major_col_ptr;
    int nb=row_major->nb;
    //Finished initialing the variables

    //Each i iteration corresponds to one row of the row major/ one col of the col major
    for(int i=0; i<nb; i++){
        //Each j iteration corresponds to one col of the row major/one row of the col major
        //The (i,j) block of the row major should be equal to the (j,i) block of the column major
        for(int j=0; j<nb; j++){

            //Accessing the row_idx and the col_ptr of the two CSCMatrixBlocked structs
            row_major_row_idx=row_major->row_idx_combined + row_major->row_idx_indices[i*nb + j];
            col_major_row_idx=col_major->row_idx_combined + col_major->row_idx_indices[j*nb + i];

            row_major_col_ptr=row_major->col_ptr_combined + row_major->col_ptr_indices[i*nb + j];
            col_major_col_ptr=col_major->col_ptr_combined + col_major->col_ptr_indices[j*nb + i];

            //Comparing the col_ptr
            for(int col_index=0; col_index<b+1; col_index++){
                //row_major_col_ptr[col_index]=-1;                  //Uncomment to deliberately create an error
                if( row_major_col_ptr[col_index] != col_major_col_ptr[col_index]){
                    printf("ERROR on compare_cols_rows on the col_ptr comparison i=%d j=%d col_index=%d, row_major_col_ptr=%d col_major_col_ptr=%d\n",i,j,col_index, row_major_col_ptr[col_index], col_major_col_ptr[col_index]);
                    exit(-1);
                }
            }

            //Comparing the row_idx
            current_non_zeros=row_major_col_ptr[b];
            for(int row_index=0; row_index<current_non_zeros; row_index++){
                //row_major_row_idx[row_index]=-1;                  //Uncomment to deliberately create an error
                if( row_major_row_idx[row_index] != col_major_row_idx[row_index]){
                    printf("ERROR on compare_cols_rows on the row_index comparison i=%d j=%d row_index=%d, row_major_row_idx=%d col_major_row_idx=%d\n",i,j,row_index, row_major_row_idx[row_index], col_major_row_idx[row_index]);
                    //exit(-1);
                }
            }



        }
    }


}


int main(int argc, char* argv[]){

    //Initialising variables:
    srand(time(NULL));
    int n=atoi(argv[1]);
    int b=atoi(argv[2]);
    if(argc<3){
        printf("Please pass the n and b parameters as arguments on blocking_test\n");
        exit(-1);
    }

    int* x=random_vector(n*n, 0.9);      //Creating some random x vector
    CSCMatrix A=array2CSC(x,n);     //Getting the CSC form 

    //Comparing the block_CSC with the fast_block_CSC to see whether the fast_block_CSC is ok
    CSCMatrix** A2=block_CSC(&A, b);
    CSCMatrixBlocked* A1=fast_block_CSC(&A, b);
    compare_block_functions(A1,A2);
    //Finished comparing the block functions


    //Comparing the two reconstruct functions, to see whether the fast_reconstruct_from_blocks is ok
    CSCMatrix* first=reconstruct_from_blocks(A2, A1->nb, n);
    CSCMatrix* fast=fast_reconstruct_from_blocks(A1);
    compare_reconstruct_functions(fast, first);
    //Finished comparing the reconstruct functions

    //Comparing the fast_block_CSC with the fast_block_CSC_cols to see whether the fast_block_CSC_cols is ok
    CSCMatrixBlocked* row_major=fast_block_CSC(&A, b);
    CSCMatrixBlocked* col_major=fast_block_CSC_cols(&A, b);
    compare_cols_rows(row_major, col_major);
    //Finished comparing the fast_block_CSC with the fast_block_CSC_cols


    free(x);
    ArrayCSCMatrixfree(A2, A1->nb);
    CSCMatrixBlocked_free(A1);
    CSCMatrixfree(first);
    CSCMatrixfree(fast);
    CSCMatrixBlocked_free(row_major);
    CSCMatrixBlocked_free(col_major);
    free(first);
    free(fast);
    free(A.col_ptr);
    free(A.row_idx);


    return 0;


}
