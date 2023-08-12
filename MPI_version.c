#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 500
//128//256 
//try on 8/16/64 cores

int main(int argc, char* argv[])
{
    int size, rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start = MPI_Wtime();

    double error, eps = 1e-5;
    double (*A)[N][N] = malloc(N * sizeof(*A));
    double (*x)[N][N] = malloc(N * sizeof(*x));
    double (*b)[N][N] = malloc(N * sizeof(*b));
    double (*x_new)[N][N] = malloc(N * sizeof(*x_new));
    int i, j, k, iter, max_iter = 1000;

    // Initialize the matrix A and vector b
    for(i = rank * N / 4; i < (rank + 1) * N / 4; i++){
        for(j = 0; j < N; j++){
            for(k = 0; k < N; k++){
                if(i == 0) {
                    x[i][j][k]=1.0; // For 3D Laplacian
                } else {
                    x[i][j][k]= 0; // boundary
                }
            }
        }
    }

    // Gauss-Seidel iteration
    for(iter = 0; iter < max_iter; iter++){
        error = 0.0;
        for(int p = 0; p < 2; p++){
            for(i = rank * N / 4; i < (rank + 1) * N / 4; i++){
                for(j = 1; j < N-1; j++){
                    for(k = (i + j + p) % 2; k < N-1; k+=2){
                        double sum = 0;
                        if(i > 0) sum += x[i-1][j][k]; 
                        if(i < N-1) sum += x[i+1][j][k];
                        if(j > 0) sum += x[i][j-1][k];
                        if(j < N-1) sum += x[i][j+1][k];
                        if(k > 0) sum += x[i][j][k-1];
                        if(k < N-1) sum += x[i][j][k+1];

                        double x_old = x[i][j][k];
                        x[i][j][k] = (b[i][j][k] + sum) / 6.0;
                        error += (x[i][j][k]-x_old) * (x[i][j][k]-x_old);
                    }
                }
            }
        }
        // Check for convergence
        double total_error;
        MPI_Allreduce(&error, &total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(total_error < eps*eps){
            if(rank == 0) {
                printf("Converged in %d iterations\n", iter);
            }
            break;
        }
        
        // Perform halo exchange
        if (rank % 2 == 0) {
            if (rank > 0) {
                MPI_Send(&x[rank * N / 4][0][0], N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&x[(rank - 1) * N / 4][0][0], N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1) {
                MPI_Send(&x[(rank + 1) * N / 4 - 1][0][0], N * N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&x[(rank + 1) * N / 4][0][0], N * N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            if (rank < size - 1) {
                MPI_Recv(&x[(rank + 1) * N / 4][0][0], N * N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&x[(rank + 1) * N / 4 - 1][0][0], N * N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
            if (rank > 0) {
                MPI_Recv(&x[(rank - 1) * N / 4][0][0], N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&x[rank * N / 4][0][0], N * N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Only rank 0 prints the solution
    if(rank == 0){
        for(i = 0; i < N; i++)
            for(j = 0; j < N; j++)
                for(k = 0; k < N; k++)
                    printf("x[%d][%d][%d] = %lf\n", i, j, k, x[i][j][k]);

        double end = MPI_Wtime();
        printf("Time taken: %f seconds\n", end - start);
    }

    free(A);
    free(x);
    free(b);
    free(x_new);

    MPI_Finalize();

    return 0;
}
