PROGRAM win_create2
    USE mpi_f08
    IMPLICIT NONE
    INTEGER::a,b,myrank
    INTEGER(KIND=MPI_ADDRESS_KIND)::size, disp
    TYPE(MPI_Win)::win
    CALL MPI_Init
    CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank)
    IF(myrank==0) a=1
    IF(myrank==1) a=10
    IF(myrank==2) a=100
    PRINT*,'myrnak:',myrank, 'a=',a
    size=STORAGE_SIZE(a)/8
    IF(MYRANK==0) THEN
        CALL MPI_Win_create(MPI_BOTTOM,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win) 
        ! CALL MPI_Win_create(a,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win)
    ELSE
        size=0
      ! 윈도우 생성 하지 않음
        ! CALL MPI_Win_create(MPI_BOTTOM,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win) 
        CALL MPI_Win_create(a,size,1,MPI_INFO_NULL,MPI_COMM_WORLD,win)
    ENDIF
    CALL MPI_Win_fence(0,win)
    IF(myrank/=0) THEN
        disp=0
        CALL MPI_GET(a,1,MPI_INTEGER,0,disp,1,MPI_INTEGER,win)
    ENDIF
    CALL MPI_Win_fence(0,win)
    IF(MYRANK==0)THEN
        disp=0
        a=200
        CALL MPI_PUT(a,1,MPI_INTEGER,1,disp,1,MPI_INTEGER,win)
        CALL MPI_PUT(a,1,MPI_INTEGER,2,disp,1,MPI_INTEGER,win)
    ENDIF
    CALL MPI_Win_fence(0,win)
    PRINT*,'MYRANK:',myrank, 'A=',a
    CALL MPI_Win_free(win)
    CALL MPI_Finalize
END PROGRAM win_create2
      
    