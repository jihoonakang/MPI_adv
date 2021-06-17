/*parallel_2D_FDM_column*/
#include <mpi.h>
#define m 6
#define n 9
void para_range(int, int, int, int, int*, int*);
int min(int, int);
main(int argc, char *argv[]){
  int i, j, nprocs, myrank ;
  double a[m][n],b[m][n];
  double works1[m],workr1[m],works2[m],workr2[m];
  int jsta, jend, jsta2, jend1, inext, iprev;
  MPI_Request isend1, isend2, irecv1, irecv2;
  MPI_Status istatus;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  para_range(0, n-1, nprocs, myrank, &jsta, &jend);
  jsta2 = jsta; jend1 = jend;
  if(myrank==0) jsta2=1;
  if(myrank==nprocs-1) jend1=n-2;
  inext = myrank + 1;
  iprev = myrank - 1;
  if (myrank == nprocs-1) inext = MPI_PROC_NULL;
  if (myrank == 0) iprev = MPI_PROC_NULL;
  for(i=0; i<m; i++)
     for(j=jsta; j<=jend; j++)  a[i][j] = i + 10.0 * j;
  if(myrank != nprocs-1) 
     for(i=0; i<m; i++) works1[i]=a[i][jend];
  if(myrank != 0) 
     for(i=0; i<m; i++) works2[i]=a[i][jsta];
  MPI_Isend(works1, m, MPI_DOUBLE, inext, 1,
            MPI_COMM_WORLD, &isend1);
  MPI_Isend(works2, m, MPI_DOUBLE, iprev, 1, 
            MPI_COMM_WORLD, &isend2);
  MPI_Irecv(workr1, m, MPI_DOUBLE, iprev, 1,
            MPI_COMM_WORLD, &irecv1);
  MPI_Irecv(workr2, m, MPI_DOUBLE, inext, 1,
            MPI_COMM_WORLD, &irecv2);
  MPI_Wait(&isend1, &istatus);
  MPI_Wait(&isend2, &istatus);
  MPI_Wait(&irecv1, &istatus);
  MPI_Wait(&irecv2, &istatus);
  if (myrank != 0)
     for(i=0; i<m; i++) a[i][jsta-1] = workr1[i];
  if (myrank != nprocs-1)
     for(i=0; i<m; i++) a[i][jend+1] = workr2[i];
  for (i=1; i<=m-2; i++)
     for(j=jsta2; j<=jend1; j++) 
         b[i][j] = a[i-1][j] + a[i][j-1] 
                      + a[i][j+1] + a[i+1][j];
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

