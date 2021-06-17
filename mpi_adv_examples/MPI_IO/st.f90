PROGRAM serial_IO2
INCLUDE 'mpif.h'
INTEGER BUFSIZE
PARAMETER (BUFSIZE = 100)
INTEGER nprocs, myrank, ierr, buf(BUFSIZE)
CHARACTER*2 number
CHARACTER*20 fname(0:16)
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
DO i = 1, BUFSIZE
    buf(i) = myrank * BUFSIZE + i
ENDDO
WRITE(number, 10) myrank
10 format(I0.2)
fname(myrank) = "testfile."//number
OPEN(UNIT=myrank+10,FILE=fname(myrank),STATUS="NEW",ACTION="WRITE")
WRITE(myrank+10,*) buf
CLOSE(myrank+10)
CALL MPI_FINALIZE(ierr)
END 

