PROGRAM cart_create
    USE mpi_f08
    IMPLICIT NONE
    INTEGER,PARAMETER::NDIM=2
    INTEGER::myrank, nprocs, dims(NDIM)=0, newprocs, newrank
    LOGICAL::reorder=.TRUE., periods(NDIM)=.false.
    TYPE(MPI_comm)::comm_cart
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD, myrank)
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs)
    CALL MPI_Dims_create(nprocs,NDIM,dims)
    if(myrank==0) PRINT*, 'DIMS=',dims
    CALL MPI_Cart_create(MPI_COMM_WORLD,NDIM,dims,periods,reorder,comm_cart)
    CALL MPI_Comm_size(comm_cart,newprocs)
    CALL MPI_Comm_rank(comm_cart,newrank);
    PRINT*,'myrank:',myrank, 'newrank:',newrank
    CALL MPI_Finalize
END PROGRAM cart_create
    