#include <stdio.h>
#include <mpi.h>

void main() {
	MPI_Init(NULL, NULL);
	printf("Hello, world!\n");
	MPI_Finalize();
}