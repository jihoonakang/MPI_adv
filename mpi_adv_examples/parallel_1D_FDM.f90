PROGRAM parallel_1D_fdm 
    INCLUDE 'mpif.h'
    PARAMETER (n=11)
    DIMENSION a(n), b(n)
    INTEGER istatus(MPI_STATUS_SIZE)
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    CALL para_range(1, n, nprocs, myrank, ista, iend)
    ista2 = ista; iend1 = iend
    IF (myrank == 0) ista2 = 2
    IF (myrank == nprocs-1) iend1 = n-1
    inext = myrank + 1; iprev = myrank-1
    IF (myrank == nprocs-1) inext = MPI_PROC_NULL
    IF (myrank == 0)        iprev = MPI_PROC_NULL
    
    DO i = ista, iend
        a(i) = 0.0d0;
        b(i) = i
    ENDDO
    
    CALL MPI_ISEND(b(iend), 1, MPI_REAL, inext, 1, &
        MPI_COMM_WORLD,  isend1, ierr)
    CALL MPI_ISEND(b(ista), 1, MPI_REAL, iprev, 1, &
        MPI_COMM_WORLD, isend2, ierr)
    CALL MPI_IRECV(b(ista-1), 1, MPI_REAL, iprev, 1, &
        MPI_COMM_WORLD, irecv1, ierr)
    CALL MPI_IRECV(b(iend+1), 1, MPI_REAL, inext, 1, &
        MPI_COMM_WORLD, irecv2, ierr)
    CALL MPI_WAIT(isend1, istatus, ierr)
    CALL MPI_WAIT(isend2, istatus, ierr)
    CALL MPI_WAIT(irecv1, istatus, ierr)
    CALL MPI_WAIT(irecv2, istatus, ierr)

    DO i = ista2, iend1
        a(i) = b(i-1) + b(i+1)
    ENDDO
    DO i = ista2, iend1
        print *, 'myrank = ',myrank,', a[',i,'] = ',a(i)
    ENDDO
    CALL MPI_FINALIZE(ierr)
END

SUBROUTINE para_range(n1, n2, nprocs, irank, ista, iend)
    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = MOD(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
    iend = ista + iwork1 - 1
    IF (iwork2 > irank) iend = iend + 1
END
  
  
   