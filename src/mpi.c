
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include "csc.h"
#include "mpi.h"

void MPI_Bcast_mat(CSCMatrixBlocked* m, int rank) {
    int nb = m->nb;
    int b = m->b;
    if (rank != 0) {
        // fill col_ptr_indices, no need for comm
        m->col_ptr_indices[0] = 0;
        for(int i=1; i<nb*nb +1; i++){
            m->col_ptr_indices[i] = (b+1)*i;
        }
    }
    
    MPI_Bcast(m->col_ptr_combined, m->nb*m->nb*(m->b+1), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m->row_idx_indices, m->nb*m->nb + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m->row_idx_combined, m->nnz , MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    CSCMatrixBlocked* A_blocked = NULL;        
    CSCMatrixBlocked* B_blocked = NULL;
    CSCMatrixBlocked* F_blocked = NULL;
    int n, nb, b, filtered;
    int a_nnz, b_nnz, f_nnz;
    CSCMatrix* C_local;

    //Initialising the environment
    MPI_Init(&argc, &argv);
    int comm_size;
    int rank;                           
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // intialization
    if (rank == 0) {
        int c;
        char* afile = NULL;
        char* bfile = NULL;
        char* ffile = NULL;

        while((c = getopt(argc, argv, "a:b:f:m:s:")) != -1) {
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
                // Block size (required)
                case 's':
                    b = atoi(optarg);
                    break;
                default:
                    fprintf(stderr, "Uknown argument, aborting");
                    abort();
            }
        }

        CSCMatrix* A = CSCfromMM(afile);
        A_blocked = block_CSC(A, b);
        CSCMatrixfree(A);
        free(A);
        CSCMatrix* B = CSCfromMM(bfile);
        B_blocked = block_CSC(B, b);
        CSCMatrixfree(B);
        free(B);
        filtered = ffile != NULL;
        if (filtered) {
            CSCMatrix* F = CSCfromMM(ffile);
            F_blocked = block_CSC(F, b);
            CSCMatrixfree(F);
            free(F);
        }

        n = A_blocked->n;
        nb = A_blocked->nb;
        a_nnz = A_blocked->nnz;
        b_nnz = A_blocked->nnz;
        if (filtered)
            f_nnz = F_blocked->nnz;
    }
    MPI_Bcast(&filtered, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nb, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&a_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (filtered)
        MPI_Bcast(&f_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("%d %d %d %d %d %d %d\n", filtered, n, nb, b, a_nnz, b_nnz, f_nnz);

    if (rank != 0) {
        A_blocked = create_empty_block_matrix(n, b, a_nnz);
        B_blocked = create_empty_block_matrix(n, b, b_nnz);
        if (filtered)
            F_blocked = create_empty_block_matrix(n, b, f_nnz);
    }
    MPI_Bcast_mat(A_blocked, rank);
    MPI_Bcast_mat(B_blocked, rank);
    if (filtered)
        MPI_Bcast_mat(F_blocked, rank);
    // end of initialization. A_blocked, B_blocked (and F_blocked) is in every process'memory now

    // bmm

    CSCMatrixBlocked_free(A_blocked);
    CSCMatrixBlocked_free(B_blocked);
    if (filtered)
        CSCMatrixBlocked_free(F_blocked);

    MPI_Finalize();
}