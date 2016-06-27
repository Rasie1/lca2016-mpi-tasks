#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MASTER 0
#define n 100000000

double dot_product(double *a, double *b, int size)
{
    double ret = 0.0;
    int i;
    for (i = 0; i < size; ++i)
        ret += a[i] * b[i];

    return ret;
}

void master(int world_rank, int world_size)
{
    // Initialize data

    double *a = (double*)malloc(n * sizeof(double));
    double *b = (double*)malloc(n * sizeof(double));

    int size = world_size - 1;
    int rank = world_rank - 1;


    int chunk_size  = (n - 1) / size + 1;
    int chunk_begin = chunk_size * rank;
    int chunk_end   = (rank + 1) * chunk_size;

    int i;
    for (i = 0; i < n; ++i)
    {
        a[i] = 1;
        b[i] = 1;
    }

    double t1, t2; 
    t1 = MPI_Wtime(); 


    // Send data for processing

    for (i = 1; i <= size; ++i)
    {
        MPI_Send(&chunk_size,
                 1,
                 MPI_INT,
                 i,
                 0,
                 MPI_COMM_WORLD);
        MPI_Send(a + (i - 1) * chunk_size, 
                 chunk_size, 
                 MPI_DOUBLE, 
                 i, 
                 0,
                 MPI_COMM_WORLD); 
        MPI_Send(b + (i - 1) * chunk_size, 
                 chunk_size, 
                 MPI_DOUBLE, 
                 i, 
                 0,
                 MPI_COMM_WORLD); 
    }

    free(a);
    free(b);


    // Gather results

    double sum = 0;

    MPI_Status  *statuses = (MPI_Status*) malloc(size * sizeof(MPI_Status));
    MPI_Request *requests = (MPI_Request*)malloc(size * sizeof(MPI_Request));
    double      *results  = (double*)     malloc(size * sizeof(double));

    for (i = 0; i < size; ++i)
    {
        int source = i + 1;
        MPI_Irecv((void*)(results + i), 
                  1, 
                  MPI_DOUBLE, 
                  source, 
                  0, 
                  MPI_COMM_WORLD, 
                  requests + i);

    }        

    MPI_Waitall(size, requests, statuses);

    for (i = 0; i < size; ++i)
    {
        sum += results[i];
    }

    printf("Dot product: %f\n", sum);
    t2 = MPI_Wtime(); 
    printf("Elapsed time is %f\n", t2 - t1); 

    free(requests);
    free(statuses);
    free(results);
}

void child(int world_rank, int world_size)
{        
    int chunk_size;
    MPI_Status status;
    MPI_Recv(&chunk_size,
             1,
             MPI_INT,
             MASTER,
             0,
             MPI_COMM_WORLD,
             &status);
    double *a = (double*)malloc(chunk_size * sizeof(double));
    double *b = (double*)malloc(chunk_size * sizeof(double));

    MPI_Recv(a,
             chunk_size,
             MPI_DOUBLE,
             MASTER,
             0,
             MPI_COMM_WORLD,
             &status);
    MPI_Recv(b,
             chunk_size,
             MPI_DOUBLE,
             MASTER,
             0,
             MPI_COMM_WORLD,
             &status);

    double result = dot_product(a, b, chunk_size);
    MPI_Send(&result, 
             1, 
             MPI_DOUBLE, 
             MASTER, 
             0,
             MPI_COMM_WORLD); 

    // printf("%f\n", result);

    free(a);
    free(b);
}

int main()
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == MASTER)
        master(world_rank, world_size);
    else
        child(world_rank, world_size);

    MPI_Finalize();

}
