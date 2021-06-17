#include <stdio.h>
#include <mpi.h>
int main(){
    int rank,nprocs,bufsize,nints;
    int buf[20]={0,},i;
    MPI_File fh;
    MPI_Offset FILESIZE;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_File_open(MPI_COMM_WORLD,"./test2.out",MPI_MODE_RDONLY,\
                    MPI_INFO_NULL,&fh);
    MPI_File_get_size(fh,&FILESIZE);
    bufsize=FILESIZE/nprocs;
    nints=bufsize/sizeof(int);
    MPI_File_read_at(fh,rank*5*sizeof(int),&buf[rank*nints],nints,\
                    MPI_INT,MPI_STATUS_IGNORE);
    printf("rank:%d  buf=",rank);
    for(i=0;i<20;i++) printf("%d  ",buf[i]);
    printf("\n");
    MPI_Finalize();
    return 0;
}
