PROGRAM get_accumulate
    USE mpi_f08
    IMPLICIT NONE
    INTEGER::a,b=0,myrank
    INTEGER(KIND=MPI_ADDRESS_KIND)::size, disp
    TYPE(MPI_Win)::win
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    IF(myrank==0) a=1
    IF(myrank==1) a=10
    IF(myrank==2) a=100
    PRINT*,'myrnak:',myrank, 'a=',a
    size=STORAGE_SIZE(a)/8
    CALL MPI_Win_create(a,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win)
    CALL MPI_Win_fence(0,win)
    
    disp=0
    IF(myrank==0) THEN
      CALL MPI_Get_accumulate(a,1,MPI_INTEGER,b,1,MPI_INTEGER,2,disp,1,&
                              MPI_INTEGER,MPI_SUM,win)
    ENDIF
    CALL MPI_Win_fence(0,win)
    PRINT*,'MYRANK:',myrank, 'A=',a
    PRINT*,'MYRANK:',myrank,'b=',b
    CALL MPI_Win_free(win)
    CALL MPI_Finalize
END PROGRAM get_accumulate
    