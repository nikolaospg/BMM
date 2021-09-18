
#include <unistd.h>

#include "csc.h"

int main(int argc, char** argv) {
    int c;
    char* afile;
    char* bfile;
    char* ffile = NULL;
    char* method;

    // http://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html
    while((c = getopt(argc, argv, "a:b:f:m:")) != -1) {
        switch (c) {
            // first bmm operand file (required)
            case 'a':
                afile = optarg;
                break;
            // second bmm operand file (required)
            case 'b':
                bfile = optarg;
                break;
            // BMM filter matrix file (optional)
            case 'f':
                ffile = optarg;
                break;
            // Method name (eg serial, block parallel, 4russians etc) (required)
            case 'm':
                method = optarg;
                break;
            default:
                fprintf(stderr, "Uknown argument, aborting");
                abort();
        }
    }

    CSCMatrix* A = CSCfromMM(afile);
    CSCMatrix* B = CSCfromMM(bfile);
    CSCMatrix* C;
    CSCMatrix* F = NULL;
    printf("f %s\n", ffile);
    if (ffile)
        F = CSCfromMM(ffile);

    C = bmm(A, B, F, method, F != NULL);

    CSCMatrixfree(A);
    CSCMatrixfree(B);
    CSCMatrixfree(C);
    if (ffile)
        CSCMatrixfree(F);
}
