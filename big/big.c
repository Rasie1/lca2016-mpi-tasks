#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#define TEST_RUNS 1
#define MASTER 0

void print_matrix(double *a, int m, int n)
{
    char* s = malloc(n * m * 3 * 2 + m + 10);
    s[0] = '\0'; 
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
            sprintf(s + strlen(s), "%.1f ", a[i * n + j]);
        sprintf(s + strlen(s), "\n");
    }
    printf("%s", s);
}

int chunk_size(int size, int parts)
{
    return (size - 1) / parts + 1;
}

int chunk_begin(int rank, int chunk_size, int parts, int procs)
{
    return chunk_size * (rank / (procs - parts + 1));
}

int chunk_end(int rank, int chunk_size, int parts, int procs)
{
    return ((rank / (procs - parts + 1)) + 1) * chunk_size;
}

int main()
{
    MPI_Init(NULL, NULL);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Generate A and B on each process

    // a: M x N, b: N x L, c: M x L
    int m = 6;
    int n = 6;
    int l = 6;

    if ((n % size) || (m % size) || (l % size))
    {
        printf("Number of processes should divide matrix sides!\n");
        return -1;
    }

    int a_w_parts = 1;
    int a_h_parts = size;
    int b_w_parts = size;
    int b_h_parts = 1;
    int c_w_parts = 1;
    int c_h_parts = size;

    int chunk_a_w = n;
    int chunk_a_h = chunk_size(m, size); 
    int chunk_b_w = chunk_size(l, size); 
    int chunk_b_h = n;
    int chunk_c_w = l;
    int chunk_c_h = chunk_size(m, size);


    int chunk_a_w_b = chunk_begin(rank, chunk_a_w, a_w_parts, size);
    int chunk_a_w_e = chunk_end  (rank, chunk_a_w, a_w_parts, size);
    int chunk_a_h_b = chunk_begin(rank, chunk_a_h, a_h_parts, size);
    int chunk_a_h_e = chunk_end  (rank, chunk_a_h, a_h_parts, size);
    int chunk_b_w_b = chunk_begin(rank, chunk_b_w, b_w_parts, size);
    int chunk_b_w_e = chunk_end  (rank, chunk_b_w, b_w_parts, size);
    int chunk_b_h_b = chunk_begin(rank, chunk_b_h, b_h_parts, size);
    int chunk_b_h_e = chunk_end  (rank, chunk_b_h, b_h_parts, size);
    int chunk_c_w_b = chunk_begin(rank, chunk_c_w, c_w_parts, size);
    int chunk_c_w_e = chunk_end  (rank, chunk_c_w, c_w_parts, size);
    int chunk_c_h_b = chunk_begin(rank, chunk_c_h, c_h_parts, size);
    int chunk_c_h_e = chunk_end  (rank, chunk_c_h, c_h_parts, size);

    double *chunk_a = (double*)malloc(sizeof(double) *
                                      chunk_a_h * chunk_a_w);
    double *chunk_b = (double*)malloc(sizeof(double) *
                                      chunk_b_h * chunk_b_w);
    double *chunk_c = (double*)malloc(sizeof(double) *
                                      chunk_c_h * chunk_c_w);

    for (int i = 0; i < chunk_a_h; ++i)
    for (int j = 0; j < chunk_a_w; ++j)
    {
        chunk_a[i * chunk_a_w + j] = 
        ((i + chunk_a_h_b) == (j + chunk_a_w_b)) ? 2.0 : 0.0;
    }
    for (int i = 0; i < chunk_b_h; ++i)
    for (int j = 0; j < chunk_b_w; ++j)
    {
        chunk_b[i * chunk_b_w + j] =// 1.0;
        ((i + chunk_b_h_b) == (j + chunk_b_w_b)) ? 4.0 : 0.0;
    }-

    char str[10000];

    double *buf = (double*)malloc(sizeof(double) * chunk_b_w * chunk_b_h);
    for (int proc = 0; proc < size; proc++)
    {
        if (rank == proc)
            for (int j = 0; j < chunk_b_h * chunk_b_w; ++j)
                buf[j] = chunk_b[j];

        MPI_Bcast(buf, 
                  chunk_b_h * chunk_b_w,
                  MPI_DOUBLE,
                  proc,
                  MPI_COMM_WORLD);

        for (int i = 0; i < chunk_c_h; ++i)
        for (int j = 0; j < chunk_c_h; ++j)
        {
            double s = 0.0;
            for (int k = 0; k < chunk_b_h; ++k)
            {
                int idx_a = k + (i) * chunk_a_w;
                int idx_b = (k) * chunk_b_w + j;
                s += chunk_a[idx_a] * buf[idx_b];
            }
            int idx_c = (i) * chunk_c_w + j + chunk_b_w * proc;
            chunk_c[idx_c] += s;
        }



    }
    printf("%i: chunk_c:\n", rank);
    print_matrix(chunk_c, chunk_c_h, chunk_c_w);
    fflush(stdout);



    // { // non-parallel
    //     int n, m, l;
    //     n = m = l = 8;
    //     double *a = (double*)malloc(sizeof(double) * n * m);
    //     double *b = (double*)malloc(sizeof(double) * n * l);
    //     double *c = (double*)malloc(sizeof(double) * m * l);


    //     for (int i = 0; i < m; ++i)
    //     for (int j = 0; j < n; ++j)
    //         a[i * n + j] = (i == j) ? 2.0 : 0.0;
    //     for (int i = 0; i < n; ++i)
    //     for (int j = 0; j < l; ++j)
    //         b[i * l + j] = (i == j) ? 4.0 : 0.0;

    //     printf("a:\n");
    //     print_matrix(a, m, n);
    //     printf("b:\n");
    //     print_matrix(b, n, l);

    //     for (int i = 0; i < m; ++i)
    //     for (int j = 0; j < l; ++j)
    //     {
    //         double s = 0.0;            
    //         for (int k = 0; k < n; ++k)
    //             s += a[n * i + k] * b[l * k + j];
    //         c[i * m + j] = s;
    //     }

    //     printf("c:\n");
    //     print_matrix(c, m, l);
    // }
    

    MPI_Finalize();
}

