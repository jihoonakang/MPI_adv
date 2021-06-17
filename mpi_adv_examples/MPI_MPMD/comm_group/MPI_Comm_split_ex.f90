PROGRAM comm_split 

    USE MPI_F08
    IMPLICIT NONE

    INTEGER(KIND=4) :: nprocs, myrank, newprocs, newrank, ikey, icolor, ranksum=0, ierr
    TYPE(MPI_COMM)  :: newcomm

    CALL MPI_INIT(ierr) 
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr) 
    
    icolor = myrank/4; ikey = mod(myrank,4) 
    
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, ikey, newcomm, ierr) 
    CALL MPI_COMM_SIZE(newcomm, newprocs, ierr) 
    CALL MPI_COMM_RANK(newcomm, newrank, ierr) 
    ! PRINT *,'myrank=', myrank, 'newprocs=', newprocs, 'newrank=',newrank 
    CALL MPI_REDUCE(newrank, ranksum, 1, MPI_INTEGER, MPI_SUM, 0, newcomm, ierr)
    PRINT *,'myrank=', myrank, 'newrank=',newrank ,'ranksum=',ranksum

    CALL MPI_FINALIZE(ierr) 
END
    