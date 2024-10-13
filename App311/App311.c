#pragma region task311
//#include <stdio.h>
//#include <mpi.h>
//
//void main() {
//	MPI_Init(NULL, NULL);
//	printf("Hello, world!\n");
//	MPI_Finalize();
//}
#pragma endregion

#pragma region task411
//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//
//void main() {
//	int ProcessRank, CommunicatorSize;
//	MPI_Init(NULL, NULL);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &CommunicatorSize);
//	printf("%d of %d\n", ProcessRank, CommunicatorSize);
//	MPI_Finalize();
//}
#pragma endregion

#pragma region task511
//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//
//void main(int argc, char* argv[]) {
//	int ProcessRank;
//	MPI_Init(argc, &argv);
//	long long N;
//	N = atoll(argv[1]);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	printf("%d:\t", ProcessRank);
//	printf("%lld\n", N);
//	MPI_Finalize();
//}
#pragma endregion

#pragma region task611

//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//
//void main(int argc, char* argv[]) {
//	int ProcessRank;
//	MPI_Init(argc, &argv);
//	long long N, i, S = 0;
//	N = atoll(argv[1]);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	for (i = 0; i < N; i++)
//		S++;
//	printf("%d:\t", ProcessRank);
//	printf("%lld\t%lld\n", N, S);
//	MPI_Finalize();
//}
#pragma endregion

#pragma region task711
//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//
//void main(int argc, char* argv[]) {
//	int ProcessRank, CommunicatorSize;
//	long long N, i, S = 0, Z;
//	MPI_Init(argc, &argv);
//	N = atoll(argv[1]);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &CommunicatorSize);
//	for (i = ProcessRank; i < N; i += CommunicatorSize)
//		S++;
//	printf("S[%d] = %lld\n", ProcessRank, S);
//	MPI_Reduce(&S, &Z, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//	if(ProcessRank == 0)
//		printf("%lld\t%lld\n", N, Z);
//	MPI_Finalize();
//}
#pragma endregion

#pragma region task811
//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//#include <time.h>
//
//void main(int argc, char* argv[]) {
//	int ProcessRank, CommunicatorSize;
//	long long N, i, S = 0, Z;
//	double T;
//	MPI_Init(argc, &argv);
//	N = atoll(argv[1]);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &CommunicatorSize);
//	clock_t t0, t1;
//	t0 = clock();
//	for (i = ProcessRank; i < N; i += CommunicatorSize)
//		S++;
//	MPI_Reduce(&S, &Z, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//	if (ProcessRank == 0) {
//		t1 = clock();
//		T = (double)(t1 - t0) / CLOCKS_PER_SEC;
//		printf("%d\t%f\t%lld\n", CommunicatorSize, T, Z);
//	}
//	MPI_Finalize();
//}

#pragma endregion

#pragma region task1011
//#include <stdio.h>
//#include <mpi.h>
//#include <stdlib.h>
//#include <time.h>
//
//void main(int argc, char* argv[]) {
//	int ProcessRank, CommunicatorSize;
//	long long N, i;
//	float T, dx, x, Z, S = 0.0;
//	MPI_Init(argc, &argv);
//	N = atoll(argv[1]);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &CommunicatorSize);
//	clock_t t0 = clock();
//	dx = 30.0 / N;
//	for (i = ProcessRank; i < N; i += CommunicatorSize) {
//		x = -10.0 + (i + 0.5) * dx;
//		S = S + 0.06 * x * x * x + 0.3 * x * x - 8.0 * x + 110.0;
//	}
//	MPI_Reduce(&S, &Z, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
//	MPI_Finalize();
//	if (ProcessRank == 0) {
//		Z = Z * dx;
//		clock_t t1 = clock();
//		T = (float)(t1 - t0) / CLOCKS_PER_SEC;
//		printf("%d\t%lld\t%.10f\t%.10f\n", CommunicatorSize, N, T, Z);
//	}
//
//}
#pragma endregion

#pragma region task1012
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

void main(int argc, char* argv[]) {
	int ProcessRank, CommunicatorSize;
	long long N, i;
	double T, dx, x, Z, S = 0.0;
	MPI_Init(argc, &argv);
	N = atoll(argv[1]);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &CommunicatorSize);
	clock_t t0 = clock();
	dx = 30.0 / N;
	for (i = ProcessRank; i < N; i += CommunicatorSize) {
		x = -10.0 + (i + 0.5) * dx;
		S = S + 0.06 * x * x * x + 0.3 * x * x - 8.0 * x + 110.0;
	}
	MPI_Reduce(&S, &Z, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Finalize();
	if (ProcessRank == 0) {
		Z = Z * dx;
		clock_t t1 = clock();
		T = (double)(t1 - t0) / CLOCKS_PER_SEC;
		printf("%d\t%lld\t%.10f\t%.10f\n", CommunicatorSize, N, T, Z);
	}

}
#pragma endregion

