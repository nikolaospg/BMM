#include "csc.h"

//Programm used to test whether the product is correct. It does this by comparing the result to the result we would have had without the CSC matrices.

int main(int argc, char* argv[]){
    srand(time(NULL));
    if(argc<3){
        printf("Give me the argument for the n and the test flag. If you want to make a filtered test, then pass 1 as the second input. If you want the non filtered test, pass 0\n");
        exit(-1);
    }
    int n=atoi(argv[1]);
    int test_flag=atoi(argv[2]);

    /*Initialising the useful matrices*/
    int* a=random_vector(n*n);          //The two operand matrices, in 1D array form
    int* b=random_vector(n*n);
    int* f=random_vector(n*n);          //the mask
    int* C_reconstructed;
    CSCMatrix A=array2CSC(a, n);         //The operand matrices, in CSC form
    CSCMatrix B=array2CSC(b, n);
    CSCMatrix F=array2CSC(f, n);
    CSCMatrix* C;

    int* c=(int*)malloc(n*n*sizeof(int));       //The result of the 1D array multiplication
    if(c==NULL){
        printf("Could not allocate memory for c, exiting\n");
        exit(-1);
    }
    /*Finished initialising*/

    /*Doing the Calculations with the CSC matrices*/
    if(test_flag==1){
        C=bmm_dsf(&A,&B,&F);                    //The result of the CSC matrix multiplication with the filter, (the function is written in csc.h)
        C_reconstructed=CSC2array(*C);        //The convertion of the CSC result to a 1D array, so we can make the test.
    }
    else{
        C=bmm_ds(&A,&B);                    //The result of the CSC matrix multiplication, (the function is written in csc.h)
        C_reconstructed=CSC2array(*C);        //The convertion of the CSC result to a 1D array, so we can make the test.
    }
    /*Finished with the CSC*/

   
    
    /*Doing the 1D array multiplication*/
    //The case of filtered BMM
    if(test_flag==1){
        int count=0;
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(f[i*n+j]!=0){
                    count=0;
                    for(int k=0; k<n; k++){
                        count=count+ a[i*n+k]*b[j+k*n];
                    }
                    c[i*n+j]=count;
                }
            }
        }
    }
    
    //The case of non filtered BMM
    else{
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
    }
    //Doing the following normalisation, so that the result is boolean
    for(int i=0; i<n*n; i++){
        if(c[i]>0){
            c[i]=1;
        }
    }
    /*Finished with the 1D multiplication*/

    /*Now the test. If you want to deliberately create an error and check whether the test is ok, uncomment the following line*/
    //c[5]=-1;    
    for(int i=0; i<n*n; i++){
        if(c[i]!=C_reconstructed[i]){
            printf("\nFor i=%d,c=%d and C=%d (C is the CSC mult. result). Error on product_test!\n",i, c[i], C_reconstructed[i]);
            exit(-1);
        }
    
    }
    /*Finished with the test*/


    free(a);
    free(b);
    free(c);
    free(f);
    free(C_reconstructed);
    CSCMatrixfree(&A);
    CSCMatrixfree(&B);
    CSCMatrixfree(&F);
    CSCMatrixfree(C);

    return 0;


}
