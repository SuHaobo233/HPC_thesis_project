#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<omp.h> 

#define N 500

int main()
{
    clock_t start, end;
    double cpu_time_used;
    
    start = omp_get_wtime();

    double error, eps = 1e-5;
    // Allocate memory on the heap
    double (*A)[N][N] = malloc(N * sizeof(*A));
    double (*x)[N][N] = malloc(N * sizeof(*x));
    double (*b)[N][N] = malloc(N * sizeof(*b));
    double (*x_new)[N][N] = malloc(N * sizeof(*x_new));
    int i, j, k, iter, max_iter = 1000;

    // Initialize the matrix A and vector b
    #pragma omp parallel for private(j,k)
    for(i = 0; i < N; i++){
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
        #pragma omp parallel for private(j,k) reduction(+:error)
        for(i = 1; i < N-1; i++){
            for(j = 1; j < N-1; j++){
                for(k = (i + j + p) % 2; k < N-1; k+=2){
                    double sum = 0;
                    if(i > 0) sum += x[i-1][j][k]; // Use the previous iteration's values
                    if(i < N-1) sum += x[i+1][j][k];
                    if(j > 0) sum += x[i][j-1][k];
                    if(j < N-1) sum += x[i][j+1][k];
                    if(k > 0) sum += x[i][j][k-1];
                    if(k < N-1) sum += x[i][j][k+1];

                    double x_old = x[i][j][k];
                    x[i][j][k] = (b[i][j][k] + sum) / 6.0; // Update in-place
                    error += (x[i][j][k]-x_old) * (x[i][j][k]-x_old);
                }
            }
        }
    }
}

        if(error < eps*eps){
            printf("Converged in %d iterations\n", iter);
            break;
        }
    }

    // Print the solution
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            for(k = 0; k < N; k++)
                printf("x[%d][%d][%d] = %lf\n", i, j, k, x[i][j][k]);

    end = omp_get_wtime();
    cpu_time_used= ((double) (end - start));

    printf("Time taken: %f seconds\n", cpu_time_used);

    // Don't forget to free the memory
    free(A);
    free(x);
    free(b);
    free(x_new);

    return 0;
}
