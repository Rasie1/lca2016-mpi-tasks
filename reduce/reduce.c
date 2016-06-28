#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MASTER 0
#define n 1000000
#define TEST_RUNS 10

double sum(double *xs, int size)
{
    double ret = 0.0;
    int i;
    for (i = 0; i < size; ++i)
        ret += xs[i];

    return ret;
}

void reduce_test(int world_rank, int world_size)
{
    double time = 0.0;
    for (int test_run = 0; test_run < TEST_RUNS; ++test_run)
    {
        double *a_master, *a_child, *results, sum;
        int chunk_size = n / world_size;
        double t1;
        if (world_rank == MASTER)
        {
            t1 = MPI_Wtime();
            a_master = (double*)malloc(n * sizeof(double));
            for (int i = 0; i < n; ++i)
                a_master[i] = 1.0;
        }
        // else
        // {
            a_child = (double*)malloc(chunk_size * sizeof(double));
        // }    
        results = (double*)malloc(chunk_size * sizeof(double));

        MPI_Scatter(a_master, 
                    chunk_size, 
                    MPI_DOUBLE,
                    a_child, 
                    chunk_size, 
                    MPI_DOUBLE, 
                    MASTER, 
                    MPI_COMM_WORLD);

        MPI_Reduce(a_child, 
                   results, 
                   chunk_size, 
                   MPI_DOUBLE,
                   MPI_SUM, 
                   MASTER, 
                   MPI_COMM_WORLD);

        if (world_rank == MASTER)
        {
            sum = 0.0;
            for (int i = 0; i < chunk_size; i++)
               sum += results[i];
            time += MPI_Wtime() - t1;
        }

    }

    if (world_rank == MASTER)
        printf("Average time for reduce is %f\n", time / (double)TEST_RUNS); 
}

void master(int world_rank, int world_size)
{
    double time = 0.0;
    int test_run, i;

    printf("Measuring simple sum...\n");

    for (test_run = 0; test_run < TEST_RUNS; ++test_run)
    {
        // Initialize data

        double *a = (double*)malloc(n * sizeof(double));

        int size = world_size - 1;
        int rank = world_rank - 1;

        int chunk_size  = (n - 1) / size + 1;
        int chunk_begin = chunk_size * rank;
        int chunk_end   = (rank + 1) * chunk_size;

        for (i = 0; i < n; ++i)
            a[i] = 1;

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
        }

        free(a);


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

        t2 = MPI_Wtime(); 

        time += t2 - t1;

        free(requests);
        free(statuses);
        free(results);
    }

    printf("Average time for simple sum is %f\n", time / (double)TEST_RUNS); 

}

void child(int world_rank, int world_size)
{        
    for (int i = 0; i < TEST_RUNS; ++i)
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

        MPI_Recv(a,
                 chunk_size,
                 MPI_DOUBLE,
                 MASTER,
                 0,
                 MPI_COMM_WORLD,
                 &status);

        double result = sum(a, chunk_size);

        MPI_Send(&result, 
                 1, 
                 MPI_DOUBLE, 
                 MASTER, 
                 0,
                 MPI_COMM_WORLD); 


        free(a);
    }
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

    reduce_test(world_rank, world_size);

    MPI_Finalize();

}
