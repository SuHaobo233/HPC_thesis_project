# Compiler
CC = gcc
MPICC = mpicc

# Flags
CFLAGS = -O2
OMPFLAGS = -fopenmp -std=c99
MPIFLAGS = -std=c99

# Target executables
TARGETS = serial openmp MPI_version

all: $(TARGETS)

serial: serial.c
	$(CC) $(CFLAGS) -o $@ $<

openmp: openmp.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $<

MPI_version: MPI_version.c
	$(MPICC) $(CFLAGS) $(MPIFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)
