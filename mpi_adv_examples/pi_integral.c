#include <stdio.h>
#include <math.h>
#include <mpi.h>
#define num_steps 1000000000

int min(int x, int y){
    int v;
    if (x>=y) v = y;
    else v = x;
    return v;
}

void para_range(int n1,int n2, int nprocs, int myrank, int *ista, int *iend){
   int iwork1, iwork2;
   iwork1 = (n2-n1+1)/nprocs; 
   iwork2 = (n2-n1+1) % nprocs; 
   *ista= myrank*iwork1 + n1 + min(myrank, iwork2);
   *iend = *ista + iwork1 - 1; 
   if(iwork2 > myrank) *iend = *iend + 1;
}

void main(int argc, char *argv[]) {
  double sum, step, x, pi;
  double tsum;
  int i, nprocs, myrank;
  int ista, iend;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  sum=0.0;
  step=1./(double)num_steps;

  para_range(1, num_steps, nprocs, myrank, &ista, &iend);
  printf(" ista=%d, iend = %d \n", ista,iend);
  for(i=ista; i<=iend; i++){
     x = (i-0.5)*step;
     sum = sum + 4.0/(1.0+x*x);
  }

  MPI_Reduce(&sum, &tsum, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

  if(myrank ==0){
    pi = step*tsum;
    printf(" numerical  pi = %.15lf \n", pi);
    printf(" analytical pi = %.15f \n", acos(-1.0));
    printf(" Error = %E \n", fabs(acos(-1.0)-pi));
  }
  MPI_Finalize();
}
