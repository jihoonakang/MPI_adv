PROGRAM Dimcreate
    USE mpi_f08
    IMPLICIT NONE
    INTEGER::dims(3)=0, myrank
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    CALL MPI_Dims_create(6,2,dims(1:2))
    if(myrank==0) PRINT*,'# of nodes (6) : dims=(',dims,')'

    dims=0
    CALL MPI_Dims_create(7,2,dims(1:2))
    if(myrank==0) PRINT*,'# of nodes (7) : dims=(',dims,')'
    dims=0

    CALL MPI_Dims_create(27,3,dims(1:3))
    if(myrank==0) PRINT*,'#of nodes (27) : dims=(',dims,')'

    CALL MPI_Finalize
END PROGRAM Dimcreate
