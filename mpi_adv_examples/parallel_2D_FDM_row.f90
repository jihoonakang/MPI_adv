PROGRAM parallel_2D_FDM_row
    INCLUDE 'mpif.h'
    PARAMETER (m = 12, n = 3)
    DIMENSION a(m,n), b(m,n)
    DIMENSION works1(n), workr1(n), works2(n), workr2(n)
    INTEGER istatus(MPI_STATUS_SIZE)
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    CALL para_range(1, m, nprocs, myrank, ista, iend)
    ista2 = ista; iend1 = iend
    IF (myrank == 0) ista2 = 2
    IF (myrank == nprocs - 1) iend1 = m-1
    inext = myrank + 1; iprev = myrank - 1
    IF (myrank == nprocs - 1) inext = MPI_PROC_NULL
    IF (myrank == 0) iprev = MPI_PROC_NULL
    
    DO j = 1, n
        DO i = ista, iend
            a(i,j) = i + 10.0 * j
        ENDDO
    ENDDO
    
    IF (myrank /= nprocs - 1) THEN
        DO j = 1, n
            works1(j) = a(iend,j)
        ENDDO
    ENDIF
    IF (myrank /= 0) THEN
        DO j = 1, n
            works2(j) = a(ista,j)
        ENDDO
    ENDIF
    CALL MPI_ISEND(works1,n,MPI_REAL,inext,1, &
        MPI_COMM_WORLD, isend1,ierr)
    CALL MPI_ISEND(works2,n,MPI_REAL,iprev,1, &
        MPI_COMM_WORLD, isend2,ierr)
    CALL MPI_IRECV(workr1,n,MPI_REAL,iprev,1, &
        MPI_COMM_WORLD, irecv1,ierr)
    CALL MPI_IRECV(workr2,n,MPI_REAL,inext,1, &
        MPI_COMM_WORLD, irecv2,ierr)
    CALL MPI_WAIT(isend1, istatus, ierr)
    CALL MPI_WAIT(isend2, istatus, ierr)
    CALL MPI_WAIT(irecv1, istatus, ierr)
    CALL MPI_WAIT(irecv2, istatus, ierr)
    IF (myrank /= 0) THEN
        DO j = 1, n
            a(ista-1,j) = workr1(j)
        ENDDO
    ENDIF

    IF (myrank /= nprocs - 1) THEN
        DO j = 1, n
            a(iend+1,j) = workr2(j)
        ENDDO
    ENDIF

    DO j = 2, n - 1
        DO i = ista2, iend1
            b(i,j) = a(i-1,j) + a(i,j-1) + a(i,j+1) + a(i+1,j)
        ENDDO
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
  