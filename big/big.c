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

int chunk_size(int size, int parts)
{
    return (size - 1) / parts + 1;
}

int chunk_begin(int rank, int chunk_size, int parts)
{
    return chunk_size * (rank / parts);
}

int chunk_end(int rank, int chunk_size, int parts)
{
    return ((rank / parts) + 1) * chunk_size;
}

int main()
{
    MPI_Init(NULL, NULL);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Generate A and B on each process

    // a: M x N, b: N x L, c: M x L
    const int m = 6;
    const int n = 6;
    const int l = 6;

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


    int chunk_a_w_b = chunk_begin(rank, chunk_a_w, a_w_parts);
    int chunk_a_w_e = chunk_end  (rank, chunk_a_w, a_w_parts);
    int chunk_a_h_b = chunk_begin(rank, chunk_a_h, a_h_parts);
    int chunk_a_h_e = chunk_end  (rank, chunk_a_h, a_h_parts);
    // int chunk_b_w_b = chunk_begin(rank, chunk_b_w, b_w_parts);
    // int chunk_b_w_e = chunk_end  (rank, chunk_b_w, b_w_parts);
    // int chunk_b_h_b = chunk_begin(rank, chunk_b_h, b_h_parts);
    // int chunk_b_h_e = chunk_end  (rank, chunk_b_h, b_h_parts);
    // int chunk_c_w_b = chunk_begin(rank, chunk_c_w, c_w_parts);
    // int chunk_c_w_e = chunk_end  (rank, chunk_c_w, c_w_parts);
    // int chunk_c_h_b = chunk_begin(rank, chunk_c_h, c_h_parts);
    // int chunk_c_h_e = chunk_end  (rank, chunk_c_h, c_h_parts);

    double *chunk_a = (double*)malloc(sizeof(double) *
                                      chunk_a_h * chunk_a_w);
    // double *chunk_b = (double*)malloc(sizeof(double) *
    //                                   chunk_b_h * chunk_b_w);
    // double *chunk_c = (double*)malloc(sizeof(double) *
    //                                   chunk_c_h * chunk_c_w);

    for (int i = chunk_a_h_b; i < chunk_a_h_e; ++i)
    for (int j = chunk_a_w_b; j < chunk_a_w_b; ++j)
        chunk_a[(i * chunk_w_b + j - ?????????????????????)] = 
        ((i + chunk_a_h_b) == (chunk_a_w_b + j)) ? 2.0 : 0.0;
    // for (int i = 0; i < n; ++i)
    // for (int j = 0; j < l; ++j)
    //     chunk_b[(i * l + j)] = 1.0;



    printf("chunk_a:\n");
    print_matrix(chunk_a, chunk_a_h, chunk_a_w);



    // non-parallel multiplication
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
    

    MPI_Finalize();
}

