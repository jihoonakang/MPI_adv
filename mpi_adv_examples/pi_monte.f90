! Example Name  : pi_monte.f90
! Compile   : $ mpif90 -g -o pi_monte.x -Wall pi_monte.f90
! Run       : $ mpirun -np 4 -hostfile hosts pi_monte.x

PROGRAM pi_monte
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER nRank, nProcs, iproc, iErr
    INTEGER :: ROOT = 0, scope = 100000000
    INTEGER :: i, nMyCount = 0, nCount = 0
    REAL :: x, y, z, pi, z1

    INTEGER status(MPI_STATUS_SIZE)

    CALL MPI_INIT(iErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, nRank, iErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, iErr)

    ! Get Random #               
    CALL INIT_RANDOM_SEED(nRank+1)
        
    DO i=0, scope-1
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(y)
        z = x*x + y*y
        z1 = SQRT(z)

        IF (z1 <= 1) nMyCount = nMyCount + 1
    END DO

    IF (nRank == ROOT) THEN
        nCount = nMyCount

        DO iproc=1, nProcs-1
            CALL MPI_RECV(nMyCount, 1, MPI_INTEGER, iproc, 55, MPI_COMM_WORLD, &
                            status, iErr)
            nCount = nCount + nMyCount
        END DO
    
        pi = 4 * REAL(nCount) / (scope * nProcs)

        WRITE (*, '(A, I2, A, I10, A)') 'Processor=', nRank, ' sending results = ', &
                nMyCount, ' to ROOT process'
        WRITE (*, '(A, I2, A, F15.10)') '# of trial is ', nProcs, ' estimate of PI is ', pi

    ELSE
        CALL MPI_SEND(nMyCount, 1, MPI_INTEGER, ROOT, 55, MPI_COMM_WORLD, iErr)
        WRITE (*, '(A, I2, A, I10)') 'Processor=', nRank, ' sending results = ', nMyCount
    END IF
    CALL MPI_FINALIZE(iErr)

    CONTAINS

    SUBROUTINE INIT_RANDOM_SEED(rank)
        IMPLICIT NONE
        INTEGER rank
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
                
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
                
        CALL SYSTEM_CLOCK(COUNT=clock)
    
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        seed = seed * rank * rank
    
        CALL RANDOM_SEED(PUT = seed)
                
        DEALLOCATE(seed)
    END SUBROUTINE
        
END
                