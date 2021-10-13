#include "bmm.h"
#include "sparse.h"


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
    //Finished initialising variables
    
    
    //Getting the blocked version and then reconstructing
    CSCMatrixBlocked* blocked= block_CSC(&A, b);
    int nb=ceil(blocked->n/((double) b));

    CSCMatrix* reconstructed=unblock_CSC(blocked, nb, A.n);
    int* y=CSC2array(*reconstructed);       //The array form of the reconstructed matrix
    //Finished reconstructing and got the final array form of the reconstructed
    
    //Testing on whether the two arrays (before and after) are the same:
    //Uncomment the following line to deliberately introduce an error
    //y[2]=-1;
    for(int i=0; i<n*n; i++){
        if(x[i]!=y[i]){
            printf("For i=%d x=%d (the initial) y=%d (the reconstructed), Error on blocking_test!\n",i, x[i], y[i]);
            exit(-1);
        }
    }
    //Finished testing


    free(y);
    CSCMatrixfree(reconstructed);
    free(reconstructed);
    free(x);
    CSCMatrixfree(&A);
    CSCMatrixBlocked_free(blocked);

    return 0;


}
