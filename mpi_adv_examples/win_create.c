#include <stdio.h>
#include <mpi.h>
int main()
{
  int a,b,myrank;
  MPI_Win win;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0) a=1;
  if(myrank==1) a=10;
  if(myrank==2) a=100;
  printf("myrank:%d   a=%d\n",myrank,a);
  MPI_Win_create(&a,sizeof(int),1,MPI_INFO_NULL,MPI_COMM_WORLD,&win);
  MPI_Win_free(&win);
  MPI_Finalize();
  return 0;
}
