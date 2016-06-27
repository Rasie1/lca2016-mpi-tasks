#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define TEST_RUNS 100
#define TEST_TIME 3
#define BUFFER_SIZE 1000000000

void process0();
void process1();

int main()
{
    MPI_Init(NULL, NULL);

    int world_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
        process0();
    else if (world_rank == 1)
        process1();

    MPI_Finalize();
}


void process0()
{
    int i;
    MPI_Status status;

    // Measure latency
    printf("Measuring latency...\n");
    double sum = 0.0;
    
    char in_buf[BUFFER_SIZE];
    char out_buf[BUFFER_SIZE];

    for (i = 0; i < TEST_RUNS; ++i)
    {
        double t1, t2; 
        t1 = MPI_Wtime(); 
        MPI_Send(out_buf,
                 0,
                 MPI_CHAR,
                 1,
                 0,
                 MPI_COMM_WORLD);
        
        MPI_Recv(in_buf,
                 0,
                 MPI_CHAR,
                 1,
                 0,
                 MPI_COMM_WORLD,
                 &status);
        t2 = MPI_Wtime(); 
        sum += t2 - t1;
    }
    printf("Latency: %f\n", sum / 2.0 / (double)TEST_RUNS); 

    // Measure average time
    printf("Measuring average time with buffer size %i...\n", BUFFER_SIZE);
    
    for (i = 0; i < TEST_RUNS; ++i)
    {
        double t1, t2; 
        t1 = MPI_Wtime(); 
        MPI_Send(out_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 1,
                 0,
                 MPI_COMM_WORLD);
        
        MPI_Recv(in_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 1,
                 0,
                 MPI_COMM_WORLD,
                 &status);
        t2 = MPI_Wtime(); 
        sum += t2 - t1;
    }
    printf("Time: %f\n", sum / 2.0 / (double)TEST_RUNS); 

    // Measure channel capacity
    printf("Measuring channel capacity...\n");
    double t = MPI_Wtime() + TEST_TIME;
    while (MPI_Wtime() <  t)
    {
        MPI_Send(out_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 1,
                 0,
                 MPI_COMM_WORLD);
    }

    // End flag
    out_buf[0] = 1;
    MPI_Send(out_buf,
             BUFFER_SIZE,
             MPI_CHAR,
             1,
             0,
             MPI_COMM_WORLD);

    printf("Capacity: %f\n", t / (double)TEST_TIME); 
}

void process1()
{
    // Measure latency
    int i;
    MPI_Status status;
    MPI_Request request;

    char in_buf[BUFFER_SIZE];
    char out_buf[BUFFER_SIZE];

    for (i = 0; i < TEST_RUNS; ++i)
    {
        MPI_Recv(in_buf,
                 0,
                 MPI_CHAR,
                 0,
                 0,
                 MPI_COMM_WORLD,
                 &status);

        MPI_Send(out_buf,
                 0,
                 MPI_CHAR,
                 0,
                 0,
                 MPI_COMM_WORLD);
    }

    for (i = 0; i < TEST_RUNS; ++i)
    {
        MPI_Recv(in_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 0,
                 0,
                 MPI_COMM_WORLD,
                 &status);

        MPI_Send(out_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 0,
                 0,
                 MPI_COMM_WORLD);
    }

    // Measure channel capacity

    while (!in_buf[0])
    {
        MPI_Recv(in_buf,
                 BUFFER_SIZE,
                 MPI_CHAR,
                 0,
                 0,
                 MPI_COMM_WORLD,
                 &status);
    }

}
