#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MASTER 0

typedef struct proc_info
{
    int rank;
    char name[MPI_MAX_PROCESSOR_NAME];
} proc_info;

void print_procs_info(int world_rank, int world_size)
{
    struct proc_info *results;
    struct proc_info result;
    results = (struct proc_info*)malloc(world_size * sizeof(struct proc_info));
    int len;
    MPI_Get_processor_name(result.name, &len);
    result.rank = world_rank;

    MPI_Datatype proc_info_type;
    MPI_Datatype types[2] = { MPI_INT, MPI_CHAR };
    int block_lengths[2] = { 1, MPI_MAX_PROCESSOR_NAME };
    MPI_Aint offsets[2];

    offsets[0] = offsetof(proc_info, rank);
    offsets[1] = offsetof(proc_info, name);

    MPI_Type_create_struct(2, 
                           block_lengths, 
                           offsets, 
                           types, 
                           &proc_info_type);
    MPI_Type_commit(&proc_info_type);

    MPI_Gather(&result, 
               1, 
               proc_info_type, 
               results, 
               1,
               proc_info_type, 
               MASTER,
               MPI_COMM_WORLD);

    if (world_rank == MASTER)
        for (int i = 1; i < world_size; ++i)
            printf("%i: \"%s\"\n", results[i].rank, results[i].name);

    MPI_Type_free(&proc_info_type);
}

int main()
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    print_procs_info(world_rank, world_size);

    MPI_Finalize();

}
