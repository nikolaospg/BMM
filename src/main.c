
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

    // timing start
    struct timespec ts_start;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    // BMM (TODO: configure MPI when needed)
    C = bmm(A, B, F, method, F != NULL);

    // timing end
    struct timespec ts_end;
    struct timespec duration;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    duration.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
    duration.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
    while (duration.tv_nsec > 1000000000) {
        duration.tv_sec++;
        duration.tv_nsec -= 1000000000;
    }
    while (duration.tv_nsec < 0) {
        duration.tv_sec--;
        duration.tv_nsec += 1000000000;
    }
    double dur_d = duration.tv_sec + duration.tv_nsec/1000000000.0;
    printf("BMM duration: %lf seconds\n", dur_d);

    CSCMatrixfree(A);
    CSCMatrixfree(B);
    CSCMatrixfree(C);
    if (ffile)
        CSCMatrixfree(F);
}
