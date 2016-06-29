#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define MASTER 0
#define MIN_WORLD_SIZE 4

void cart_test(MPI_Comm comm)
{
    int rank;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm cart_comm;

    printf("Cart shift: comm size %i, proc rank %i\n", size, rank);

    int cyclic = 1;

    MPI_Cart_create(comm, 
                    1, 
                    &size, 
                    &cyclic, 
                    1,
                    &cart_comm);

    int source, dest;
    MPI_Cart_shift(cart_comm, 0, 1, &source, &dest);

    int send_value = rank + 100;
    int recv_value;
    MPI_Status status;
    MPI_Sendrecv(&send_value,
                 1,
                 MPI_INT,
                 dest,
                 0,
                 &recv_value,
                 1,
                 MPI_INT,
                 source,
                 0,
                 cart_comm,
                 &status);
    printf("C%i: sent %i to C%i, received %i from C%i\n", 
           rank, send_value, dest, recv_value, source);

}

void graph_test(MPI_Comm comm)
{
    int rank;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm graph_comm;

    printf("Graphs: comm size %i, proc rank %i\n", size, rank);

    int *neighbors_count = (int*)malloc(sizeof(int) * size);
    int *edges = (int*)malloc(sizeof(int) * 2 * size);
    
    neighbors_count[MASTER] = size - 1;
    for (int i = 1; i < size; ++i)
        neighbors_count[i] = neighbors_count[i - 1] + 1;

    int j = 0;
    for (int i = 1; i < size; ++i)
        edges[j++] = i;
    for (int i = 1; i < size; ++i)
        edges[j++] = MASTER;

    MPI_Graph_create(comm, 
                     size, 
                     neighbors_count, 
                     edges, 
                     1, 
                     &graph_comm);

    int *neighbors = (int*)malloc(size * sizeof(int)); 
    MPI_Graph_neighbors(graph_comm, rank, size - 1, neighbors);
    int current_neighbor_count;
    MPI_Graph_neighbors_count(graph_comm, rank, &current_neighbor_count);


    int send_value = rank + 100;
    int recv_value;
    MPI_Status status;

    for (int i = 0; i < current_neighbor_count; ++i)
    {
        MPI_Sendrecv(&send_value,
                     1,
                     MPI_INT,
                     neighbors[i],
                     0,
                     &recv_value,
                     1,
                     MPI_INT,
                     neighbors[i],
                     0,
                     graph_comm,
                     &status);
        printf("g%i: sent %i to g%i, received %i from g%i\n", 
               rank, send_value, neighbors[i], recv_value, neighbors[i]);
    }

    free(neighbors);
    free(neighbors_count);
    free(edges);
}

void topologies(int world_rank, int world_size)
{
    MPI_Comm comm;
    int cart_process = world_rank % 2;

    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, cart_process, world_rank, &comm);

    if (cart_process)
        cart_test(comm);
    else
        graph_test(comm);

}

int main()
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size < MIN_WORLD_SIZE)
    {
        printf("Impossible to run with fewer than 4 processes");
        return -1;
    }

    topologies(world_rank, world_size);



    MPI_Finalize();

}
