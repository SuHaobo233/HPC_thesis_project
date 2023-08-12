# MATRIX MULTIPLICATION IN DISTRIBUTED SYSTEMS
This project encompasses several parallel and serial versions of programs, implemented using MPI, OpenMP, and standard serial programming techniques.

**File List**

**1**. Dockerfiles
The project contains several Dockerfiles meant for setting up execution environments.

- `Dockerfile.mpi`: Establishes an environment that supports MPI (Message Passing Interface).
- `Dockerfile.openmp`: Sets up an environment that supports OpenMP.
- `Dockerfile.serialcode`: Creates an environment for executing serial code.

**2**. Code Files
- `MPI_version.c`: A parallel version of the code implemented using MPI.
- `openmp.c`: A parallel version of the code implemented using OpenMP.
- `serial.c`: A serial version of the code.

**3**. Makefile
A generic `Makefile` for compiling the aforementioned C code files.

## Compilation and Execution

### Compilation
Compile all the source code files using:
```
$make
```

### Running with Docker
To use the Dockerfiles and run the programs:

1. Build the Docker image:

For MPI:
```
$docker build -t mpi_image -f Dockerfile.mpi .
```
For Opnemp:
```
$docker build -t openmp_image -f Dockerfile.openmp .
```
For serial code:
```
$docker build -t serial_image -f Dockerfile.serialcode .
```
2.Run the Docker container:

For MPI:
```
$docker run mpi_image mpirun -np 4 ./MPI_version
```
For Opnmp:
```
$docker run openmp_image ./openmp
```

For serial code:
```
$docker run serial_image ./serial
```

