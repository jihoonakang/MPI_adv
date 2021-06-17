#include <stdio.h>
#include "mpi.h"
int main()
{
    MPI_File fh;
    int buf[20]={0,}, rank,i,bufsize,nints,offset,data[20]={0,};
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    for(i=0;i<20;i++) buf[i]=i+1;
    MPI_File_open(MPI_COMM_WORLD,"test2.out",\
                    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
    offset=rank*(20/4)*sizeof(int); // 20: # of data, 4: # of processes
    MPI_File_seek(fh,offset,MPI_SEEK_SET);
    MPI_File_write(fh,&buf[rank*5],5,MPI_INT,MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
    MPI_Finalize();
    return 0;
}
