
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include "bmm.h"
#include "sparse.h"
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
    int test; // used for validation
    CSCMatrix* C_actual = NULL; //used for validation
    struct timespec ts_start, ts_end;
    int n, nb, b, filtered;
    int a_nnz, b_nnz, f_nnz;

    int comm_size;
    int rank; 
    int nnz_tag = 1;
    int row_idx_tag = 2;
    int col_ptr_tag = 3;
    MPI_Status status;

    //Initialising the environment
    MPI_Init(&argc, &argv);                          
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // intialization
    if (rank == 0) {
        int c;
        char* afile = NULL;
        char* bfile = NULL;
        char* ffile = NULL;
        test = 0;

        while((c = getopt(argc, argv, "a:b:f:m:s:t")) != -1) {
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
                    nb = atoi(optarg);
                    break;
                // Whether to run test after result (optional)
                case 't':
                    test = 1;
                    break;
                default:
                    fprintf(stderr, "Uknown argument, aborting");
                    abort();
            }
        }

        CSCMatrix* A = CSCfromMM(afile);
        CSCMatrix* B = CSCfromMM(bfile);
        filtered = ffile != NULL;
        CSCMatrix* F = NULL;
        if (filtered) {
            F = CSCfromMM(ffile);
        }
        if (test) {
            if (filtered)
                C_actual = bmm_ssf(A, B, F);
            else
                C_actual = bmm_ss(A, B);
        }
        printf("%s, nb=%d, filtered=%d, comm_size=%d\n", afile, nb, filtered, comm_size);
        clock_gettime(CLOCK_MONOTONIC, &ts_start); // time start

        b = ceil(A->n/((double)nb));
        A_blocked = block_CSC(A, b);
        CSCMatrixfree(A);
        free(A);
        B_blocked = block_CSC(B, b);
        CSCMatrixfree(B);
        free(B);
        if (filtered) {
            F_blocked = block_CSC(F, b);
            CSCMatrixfree(F);
            free(F);
        }
        struct timespec ts_end_blocking;
        clock_gettime(CLOCK_MONOTONIC, &ts_end_blocking);
        double dur_d = duration_secs(ts_start, ts_end_blocking);
        printf("block duration: %lf seconds\n", dur_d);

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
    // end of initialization. A_blocked, B_blocked (and F_blocked) is in every process' memory now

    // bmm
    // each process will be dealt a block of C, in row major fashion
    // each process will process ceil(nb^2/comm_size) blocks at index rank+i*comm_size, i=0,..,ceil(nb^2/comm_size)-1
    int max_blocks_per_process = ceil(nb*nb/((double)comm_size));
    CSCMatrix** process_blocks = malloc(max_blocks_per_process*sizeof(CSCMatrix*));
    #pragma omp parallel for if(!filtered) // non filtered BMM is serial, use OpenMP across block BMMs
    for (int i = 0; i < max_blocks_per_process; i++) {
        int block_index = rank+i*comm_size;
        if (block_index >= nb*nb) continue; // no more blocks
        int p = block_index/nb;
        int q = block_index%nb;

        // Cp,q = 0
        process_blocks[i] = malloc(sizeof(CSCMatrix));
        process_blocks[i]->n = b;
        process_blocks[i]->col_ptr = calloc(b+1, sizeof(int));
        process_blocks[i]->row_idx = NULL; // NULL can be passed to realloc later
        if (filtered)
            block_bmmf(A_blocked, B_blocked, F_blocked, &process_blocks[i], p, q);
        else
            block_bmm(A_blocked, B_blocked, &process_blocks[i], p, q);
    }
    // bmm on each block computed
    // free A_blocked, B_blocked, no longer needed
    CSCMatrixBlocked_free(A_blocked);
    CSCMatrixBlocked_free(B_blocked);
    if (filtered)
        CSCMatrixBlocked_free(F_blocked);

    // send blocks to coordinator and reconstruct C
    if (rank == 0) {
        CSCMatrix** all_blocks = malloc(nb*nb*sizeof(CSCMatrix*));
        int* nnz_map = malloc(nb*nb*sizeof(int)); // map that maps block to NNZ
        // fill nnz_map
        for (int i = 0; i < nb*nb; i++) {
            int block_rank = i%comm_size; // rank that computed block
            if (block_rank == 0) { // get from self
                int block_idx = i/comm_size;
                CSCMatrix* block = process_blocks[block_idx];
                nnz_map[i] = block->col_ptr[block->n];
            } else { //receive from other process
                MPI_Recv(&nnz_map[i], 1, MPI_INT, block_rank, nnz_tag, MPI_COMM_WORLD, &status);
            }
        }
        // nnz_map filled, now get matrices
        for (int i = 0; i < nb*nb; i++) {
            int block_rank = i%comm_size; // rank that computed block
            all_blocks[i] = malloc(sizeof(CSCMatrix));
            all_blocks[i]->n = b;
            all_blocks[i]->col_ptr = malloc((b+1)*sizeof(int));
            all_blocks[i]->row_idx = malloc(nnz_map[i]*sizeof(int));
            if (block_rank == 0) { // get from self
                int block_idx = i/comm_size;
                CSCMatrix* block_src = process_blocks[block_idx];
                memcpy(all_blocks[i]->col_ptr, block_src->col_ptr, (b+1)*sizeof(int));
                memcpy(all_blocks[i]->row_idx, block_src->row_idx, nnz_map[i]*sizeof(int));
            } else { //receive from other process
                MPI_Recv(all_blocks[i]->col_ptr, b+1, MPI_INT, block_rank, col_ptr_tag, MPI_COMM_WORLD, &status);
                MPI_Recv(all_blocks[i]->row_idx, nnz_map[i], MPI_INT, block_rank, row_idx_tag, MPI_COMM_WORLD, &status);
            }
        }
        // block matrix attained, now convert back to CSC
        struct timespec ts_start_unblock;
        clock_gettime(CLOCK_MONOTONIC, &ts_start_unblock);
        int total_nnz = 0;
        for (int i = 0; i<nb*nb; i++) total_nnz+=nnz_map[i];
        CSCMatrixBlocked* C_blocked = blockcsc_tobsc(all_blocks, nb, n, total_nnz);

        free(nnz_map);
        for (int i = 0; i < nb*nb; i++) {
            CSCMatrixfree(all_blocks[i]);
            free(all_blocks[i]);
        }
        free(all_blocks);

        CSCMatrix* C = unblock_CSC(C_blocked, nb, n);
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        double dur_d = duration_secs(ts_start_unblock, ts_end);
        printf("unblock duration: %lf seconds\n", dur_d);
        dur_d = duration_secs(ts_start, ts_end);
        printf("MPI BMM duration: %lf seconds\n", dur_d);
        if (test) {
            csc_validate(C, C_actual);
            CSCMatrixfree(C_actual);
            free(C_actual);
        }
    } else {
        // send NNZ
        for (int i = 0; i < max_blocks_per_process; i++) {
            int block_index = rank+i*comm_size;
            if (block_index >= nb*nb) break; // no more blocks
            CSCMatrix* block = process_blocks[i];
            MPI_Send(&block->col_ptr[block->n], 1, MPI_INT, 0, nnz_tag, MPI_COMM_WORLD);
        }
        // send blocks
        for (int i = 0; i < max_blocks_per_process; i++) {
            int block_index = rank+i*comm_size;
            if (block_index >= nb*nb) break; // no more blocks
            CSCMatrix* block = process_blocks[i];
            MPI_Send(block->col_ptr, block->n+1, MPI_INT, 0, col_ptr_tag, MPI_COMM_WORLD);
            MPI_Send(block->row_idx, block->col_ptr[block->n], MPI_INT, 0, row_idx_tag, MPI_COMM_WORLD);
        }
    }
    // free other common resources
    for (int i = 0; i < max_blocks_per_process; i++) {
        int block_index = rank+i*comm_size;
        if (block_index >= nb*nb) break; // no more blocks
        CSCMatrixfree(process_blocks[i]);
        free(process_blocks[i]);
    }
    free(process_blocks);

    MPI_Finalize();
}
