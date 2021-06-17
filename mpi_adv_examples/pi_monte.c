/*
	Example Name	: pi_monte.c
	Compile	: $ mpicc -g -o pi_monte -Wall pi_monte.c
	Run		: $ mpirun -np 4 -hostfile hosts pi_monte
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define	SCOPE	100000000

/*
** PI (3.14) Monte Carlo Method
*/
int main(int argc, char *argv[])
{
	int nProcs, nRank, proc, ROOT = 0;
	MPI_Status status;
	int nTag = 55;

	int i, nCount = 0, nMyCount = 0;
	double x, y, z, pi, z1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	srand(time(NULL)*(nRank+1));

	for (i = 0; i < SCOPE; i++) {
		x = (double)rand() / (double)(RAND_MAX);
		y = (double)rand() / (double)(RAND_MAX);

		z = x*x + y*y;
		z1 = sqrt(z);

		if (z1 <= 1)	nMyCount++;
	}

	if (nRank == ROOT) {
		nCount = nMyCount;

		for (proc=1; proc<nProcs; proc++) {
			MPI_Recv(&nMyCount, 1, MPI_REAL, proc, nTag, MPI_COMM_WORLD,
                             &status);
			nCount += nMyCount;
		}

		pi = (double)4*nCount/(SCOPE * nProcs);
		printf("Processor %d sending results = %d to ROOT processor\n", nRank,
                 nMyCount);
		printf("\n # of trials (cpu#: %d, time#: %d) = %d, estimate of pi is %f\n",
				nProcs, SCOPE, SCOPE*nProcs, pi);

	}
	else {
		printf("Processor %d sending results = %d to ROOT processor\n", nRank,
                 nMyCount);
		MPI_Send(&nMyCount, 1, MPI_REAL, ROOT, nTag, MPI_COMM_WORLD);
	}

	printf("\n");

	MPI_Finalize();

	return 0;
}
