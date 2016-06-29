#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define TEST_RUNS 1
#define MASTER 0

void print_matrix(double *a, int m, int n)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
            printf("%.1f ", a[i * n + j]);
        printf("\n");
    }
}

int main()
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    if (world_rank == MASTER)
    {
        const int m = 6;
        const int n = 6;
        const int l = 6;

        double *a = (double*)malloc(sizeof(double) * n * m);
        double *b = (double*)malloc(sizeof(double) * n * l);
        double *c = (double*)malloc(sizeof(double) * m * l);


        for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            a[i * n + j] = (i == j) ? 2.0 : 0.0;
        for (int i = 0; i < n; ++i)
        for (int j = 0; j < l; ++j)
            b[i * l + j] = 1.0;

        printf("a:\n");
        print_matrix(a, m, n);
        printf("b:\n");
        print_matrix(b, n, l);

        for (int i = 0; i < m; ++i)
        for (int j = 0; j < l; ++j)
        {
            double s = 0.0;            
            for (int k = 0; k < n; ++k)
                s += a[n * i + k] * b[l * k + j];
            c[i * m + j] = s;
        }

        printf("c:\n");
        print_matrix(c, m, l);
    }

    MPI_Finalize();
}

