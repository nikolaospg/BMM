CC=gcc
MPICC=mpicc
CPPC=g++

BINS_DIR=bin
MMIO_LIB=$(BINS_DIR)/mmio.o
C_SOURCES = mmio/mmio.c
SOURCES = $(C_SOURCES) $(CPP_SOURCES)

CFLAGS=-Wall -O3 -Iinclude -Immio -std=c99 -D_POSIX_C_SOURCE=199309L -fopenmp # https://stackoverflow.com/questions/29666937/error-clock-monotonic-undeclared-first-use-in-this-function
DEBUG_CFLAGS=-Wall -g -fsanitize=address -Iinclude -Immio -std=c99 -D_POSIX_C_SOURCE=199309L -fopenmp
CPPFLAGS=-std=c++0x $(CFLAGS)
OMPFLAGS=-std=c++0x -Wall -O3 -fopenmp -DOMP
PTHREADSFLAGS=-pthread -DPTHREADS
LDFLAGS=$(MMIO_LIB) -lm

default: all

.PHONY: clean

bin:
	mkdir -p $@

data:
	chmod +x get_data.bash
	./get_data.bash

mmio: bin
	$(CC) -c $(CFLAGS) -o $(MMIO_LIB) $(C_SOURCES)

bmm: bin | mmio
	$(CC) $(CFLAGS) -o $(BINS_DIR)/$@ src/bmm.c $(LDFLAGS)

mpi: bin | mmio
	$(MPICC) $(CFLAGS) -o $(BINS_DIR)/$@ src/mpi.c $(LDFLAGS)

test: bin | mmio mpi
	$(CC) $(CFLAGS) -o $(BINS_DIR)/product_test test/product_test.c $(LDFLAGS)
	$(CC) $(CFLAGS) -o $(BINS_DIR)/conversion_test test/conversion_test.c $(LDFLAGS)
	$(CC) $(CFLAGS) -o $(BINS_DIR)/blocking_test test/blocking_test.c $(LDFLAGS) 
	#$(CC) $(CFLAGS) -o $(BINS_DIR)/demo_bss test/demo_bss.c $(LDFLAGS) 
	chmod +x all_tests.bash
	./all_tests.bash

bash: bin
	chmod +x bmm_self.bash
	chmod +x bmm_self_mpi.bash
	chmod +x benchmark.bash
	chmod +x benchmark_mpi.bash

all: bmm test bash mpi

clean:
	rm -rf $(BINS_DIR)
