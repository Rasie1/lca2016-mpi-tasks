#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MASTER 0

struct proc_info
{
    int rank;
    const char* name;
};

void print_procs_info(int world_rank, int world_size)
{
    proc_info *results;
    if (world_rank == MASTER)
    {
        results = (proc_info*)malloc(world_size * sizeof(proc_info));
    }

    

}

int main()
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    reduce_test(world_rank, world_size);

    MPI_Finalize();

}
