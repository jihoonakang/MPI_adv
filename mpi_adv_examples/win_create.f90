PROGRAM win_create
    USE mpi_f08
    IMPLICIT NONE
    INTEGER::a,b,myrank
    INTEGER(KIND=MPI_ADDRESS_KIND)::size,lb
    TYPE(MPI_Win)::win
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    IF(myrank==0) a=1
    IF(myrank==1) a=10
    IF(myrank==2) a=100
    PRINT*,'myrnak:',myrank, 'a=',a
    size=STORAGE_SIZE(a)/8
    CALL MPI_Type_get_extent(MPI_INTEGER,lb,size)
    CALL MPI_Win_create(a,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win)
    CALL MPI_Win_free(win)
    CALL MPI_Finalize
END PROGRAM win_create
    