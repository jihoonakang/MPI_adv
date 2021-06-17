PROGRAM cart_shift
    USE mpi_f08
    IMPLICIT NONE
    TYPE(MPI_Comm)::oldcomm, newcomm
    INTEGER::ndims=2, dimsize(0:1)
    LOGICAL::periods(0:1), reorder
    INTEGER::myrank,nprocs,i,j,rank
    INTEGER::coords(0:1)
    INTEGER::direction, disp,src,dest
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    CALL MPI_Comm_size(MPI_COMM_WORLD,nprocs)
    oldcomm=MPI_COMM_WORLD
    dimsize=(/3,2/)
    periods=(/.TRUE., .FALSE./)
    reorder=.TRUE.
    CALL MPI_Cart_create(oldcomm, ndims,dimsize,periods, reorder,newcomm)
    direction=0;  disp=1
    CALL MPI_Cart_shift(newcomm, direction,disp,src,dest)
    PRINT*,'rank:',myrank,'source=',src,'destination=',dest
    CALL MPI_Finalize
END PROGRAM cart_shift
    