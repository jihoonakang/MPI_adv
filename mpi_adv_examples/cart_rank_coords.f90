PROGRAM cart_rank_coords
    USE mpi_f08
    IMPLICIT NONE
    TYPE(MPI_Comm)::oldcomm, newcomm
    INTEGER::ndims=2, dimsize(0:1)
    LOGICAL::periods(0:1), reorder
    INTEGER::myrank,nprocs,i,j,rank
    INTEGER::coords(0:1)
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    CALL MPI_Comm_size(MPI_COMM_WORLD,nprocs)
    oldcomm=MPI_COMM_WORLD
    dimsize=(/3,2/)
    periods=(/.TRUE., .FALSE./)
    reorder=.TRUE.
    CALL MPI_Cart_create(oldcomm, ndims,dimsize,periods, reorder,newcomm)
    IF(myrank==0)THEN
      DO i=0,dimsize(0)-1
        DO j=0,dimsize(1)-1
          coords=(/i,j/)
          CALL MPI_Cart_rank(newcomm,coords,rank)
          PRINT*,'coords=',coords,'rank=',rank
        END DO
      END DO
    ENDIF
    IF(myrank==0)THEN
      DO rank=0,nprocs-1
        CALL MPI_Cart_coords(newcomm,rank,ndims,coords);
        PRINT*,'rank=',rank,'coords=',coords
      END DO
    END IF
    CALL MPI_Finalize
END PROGRAM cart_rank_coords
      