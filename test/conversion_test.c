/*This is a script demonstrating that the convertion from a row major to a CSC and back is correct. Just run the script.sh and you will see.*/


#include "algorithms.h"


int main(int argc, char* argv[]){

    srand(time(NULL));
    if(argc<2){                                         //The n number is given as a command line argument
        printf("Give me the n argument\n");
        exit(-1);
    }
    int n=atoi(argv[1]);

    /*Creating the x array, getting the equivalent CSC struct, and then converting this as well*/
    int* x=random_vector(n*n);      //The initial array
    CSCArray shit=array2CSC(x, n);   //The CSC struct we get 
    int * y= CSC_convert(shit);     //The final array we get, after converting the converted one 
    /*Finished with the convertions*/

    /*Doing the test*/
    //x[5]=3;       //Uncomment to deliberately create an error and see that the check is ok
    for(int i=0; i<n*n; i++){
        if(x[i]!=y[i]){
            printf("For i=%d, x[i]=%d, y[i]=%d\n. Error!",i,x[i],y[i]);
            exit(-1);
        }
    
    }
    /*Finished with the test*/
    CSCArrayfree(&shit);
    free(y);
    free(x);
    return 0;
}

