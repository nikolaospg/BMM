# BMM
## Building instructions
1. Download the data from SuiteSparse matrix collection  
`make data`  
You can change `dataset.txt` if you want to download a different set of matrices. However you **will** need mycielskian13 to be able to run `all_tests.bash`
2. Build `bmm` and `mpi` targets  
`make bmm`  
`make mpi`
## Running
We have two main executables, `bmm` and `mpi` both of which are generated into the `bin` folder. Use `mpi` to run the MPI implementation of BMM and `bmm`
to try the various different BMM algorithms.
## Usage

```
./bin/bmm -a <matrix A> -b <matrix B> (-f <matrix F>) -m <method> (-s <numblocks>)
```

```
mpiexec -np <numprocesses> ./bin/mpi -a <matrix A> -b <matrix B> (-f <matrix F>) (-s <numblocks>) (-t)
```

* Matrices `A`, `B` and `F` are the operands of the BMM operation. Matrix F is optional.
* `numblocks` is the number of blocks across each matrix dimension. The total number of blocks will be `numblocks*numblocks`
* `-t` is an optional flag for `mpi` that validates the output by comparing it to a serial BMM implementation
* `method` is the block-level BMM implementation that is going to be used. It can be one of the methods described in the next section
### BMM methods
* `ds`  
Serial implementation using the binary matrix multiplication definition, very slow.
* `cp`  
CSR*CSC style implementation parallelized in OpenMP. Fast for filtered BMM, but for low number of threads and very sparse matrices such as `europe_osm`,
filtered `ss` is likely faster. Not supported for non-filtered BMM, use `ss` instead.
* `ss`  
symbolic, SMMP-like BMM with only one pass. Very fast for non-filtered BMM, but mostly outclassed by `cp` for filtered BMM, especially on denser matrices.
* `sp`  
An attempt to parallelize `ss`, which ended up needing too many threads to become faster than `ss`.
* `a`  
Adaptive method which chooses between `ss` and `cp` based on matrix desnity. Only supported for filtered BMM.
* `bp`  
Block BMM with OpenMP parallelization. Uses `ss` for non-filtered BMM and `a` for filtered BMM.
### Convenient bash scripts
```
./bmm_self.bash <matrix> <method> (-f)
```
Performs BMM with A=B=F=<matrix> so you don't have to type the same matrix 3 times. Cannot be used with `bp`
```
./bmm_self_mpi.bash <matrix> <blocksize> (-f)
```
Similar to the previous one but for `mpi`. Instead of `method` you specify `blocksize`. Default number of processes is 2, if you need more, you can edit the file directly
```
./benchmark.bash <method> (-f)
```
Run `bmm_self.bash` for multiple matrices in a row, specified in `benchmark.bash`
  ```
./benchmark_mpi.bash <blocksize> (-f)
```
Similar to the previous one but for `mpi`. Instead of `method` you specify `blocksize`
## Run tests
All tests are in the `test` folder, which cover every BMM implementation. To run all tests, just run
```
make test
```
