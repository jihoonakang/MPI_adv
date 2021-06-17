PROGRAM neighborhood_allgather
    USE mpi_f08
    IMPLICIT NONE
    INTEGER::myrank,procs
    TYPE(MPI_Comm)::comm_cart
    INTEGER::NDIM=2,dims(2), coords(2),recvbuf(4)=-10
    LOGICAL::periods(2)=.FALSE., reorder=.FALSE.
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    CALL MPI_Comm_size(MPI_COMM_WORLD,procs)
    dims=(/4,4/)
    CALL MPI_Cart_create(MPI_COMM_WORLD,NDIM,dims,periods,reorder,comm_cart)
    CALL MPI_Cart_coords(comm_cart,myrank,NDIM,coords)
    CALL MPI_Neighbor_allgather(myrank,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,comm_cart)
    PRINT*,'myrank:',myrank,'recvbuf=',recvbuf
    CALL MPI_Finalize
END PROGRAM neighborhood_allgather
    