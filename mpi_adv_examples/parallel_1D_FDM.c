/*parallel_1D_fdm*/
#include <mpi.h>
#define n 11
void para_range(int, int, int, int, int*, int*);
int min(int, int);
main(int argc, char *argv[]){
    int i, nprocs, myrank ;
    double a[n], b[n];
    int ista, iend, ista2, iend1, inext, iprev;
    MPI_Request isend1, isend2, irecv1, irecv2;
    MPI_Status istatus;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    para_range(0, n-1, nprocs, myrank, &ista, &iend);
    ista2 = ista; iend1 = iend;
    if(myrank==0) ista2=1;
    if(myrank==nprocs-1) iend1=n-2;
    inext=myrank+1; iprev=myrank-1;
    if(myrank==nprocs-1) inext=MPI_PROC_NULL;
    if(myrank==0) iprev=MPI_PROC_NULL; 

    for(i=ista; i<=iend; i++) {
        b[i] = i;
        a[i] = 0.0;
    }
    MPI_Isend(&b[iend], 1, MPI_DOUBLE, inext ,1, MPI_COMM_WORLD,  &isend1);
    MPI_Isend(&b[ista], 1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &isend2);
    MPI_Irecv(&b[ista-1], 1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &irecv1);
    MPI_Irecv(&b[iend+1], 1, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &irecv2);
    MPI_Wait(&isend1, &istatus);
    MPI_Wait(&isend2, &istatus);
    MPI_Wait(&irecv1, &istatus);
    MPI_Wait(&irecv2, &istatus);
    for(i=ista2; i<=iend1; i++) a[i] = b[i-1] + b[i+1];

    for(i=ista2; i<=iend1; i++) {
        printf("myrank = %d, a[%d] = %f\n",myrank, i, a[i]);
    }
    MPI_Finalize();
}

void para_range(int n1,int n2, int nprocs, int myrank, int *ista, int *iend){
   int iwork1, iwork2;
   iwork1 = (n2-n1+1)/nprocs; 
   iwork2 = (n2-n1+1) % nprocs; 
   *ista= myrank*iwork1 + n1 + min(myrank, iwork2);
   *iend = *ista + iwork1 - 1; 
   if(iwork2 > myrank) *iend = *iend + 1;
}

int min(int x, int y){
    int v;
    if (x>=y) v = y;
    else v = x;
    return v;
}

