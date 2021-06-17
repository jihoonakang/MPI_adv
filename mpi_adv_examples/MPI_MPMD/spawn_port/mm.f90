! reference : https://www.hpc.cineca.it/content/mm-fortran-mpi

program mm

    use mpi
    implicit none
    integer :: n, size
    integer :: nprocs, ierr, myrank, status(MPI_STATUS_SIZE)
    real(kind(1.d0)), dimension(:,:), allocatable   :: a, b, c
    real(kind(1.d0)) :: d
    real(kind(1.d0)) :: time1, time2, gflops, pi

    integer            :: i, j, k, jstart, jend
    character(len=128) :: command
    character(len=80) :: arg
    character(len=100) :: portname

    integer :: MPI_COMM_MASTER, MPI_COMM_WORKER0

!mpi stuff
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    call MPI_Comm_get_parent(MPI_COMM_MASTER,ierr)

!reading matrix size
    if(myrank.eq.0) then
        write(*,*) "[mm] MPI tasks = ", nprocs
        call get_command_argument(0,command)
        if (command_argument_count() /= 1) then
            write(0,*) '[mm] Usage:', trim(command), '   matrix size'
            call MPI_ABORT(MPI_COMM_WORLD, 911,ierr)
            call MPI_FINALIZE(ierr)
            stop
        else
            call get_command_argument(1,arg)
            read(arg,*) n
        endif
    endif
    write(*,*) "[mm] MPI ranks = ", myrank
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!a check
    if (n > 0 ) then
        if (myrank == 0 ) then
            write(*,*) '[mm] Matrix size is ', n
        endif
    else
        write(*,*) "[mm] Error, matrix size is ", n
        call MPI_FINALIZE(ierr)
        stop
    endif
!
    if(nprocs /= 2) then
        if(myrank == 0) then
            write(*,*) "[mm] Error, nprocs=", nprocs, "is not 2!!!!"
            write(*,*) "[mm] This (stupid) code works only with 2 task!!!"
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
    endif
!
!allocation/inizializations
    allocate(a(n,n),stat=ierr)
    if(ierr/=0) STOP '[mm] a matrix allocation failed'
    allocate(b(n,n),stat=ierr)
    if(ierr/=0) STOP '[mm] b matrix allocation failed'
    allocate(c(n,n),stat=ierr)
    if(ierr/=0) STOP '[mm] c matrix allocation failed'
    c = 0.d0

    if(myrank == 0) then
        call random_number(a)
        call random_number(b)
    endif
!
    time1 = MPI_Wtime()
!sending a and b elements
    if(myrank.eq.1) then
        call mpi_recv(a(1,1), n*n, MPI_DOUBLE,0,1, MPI_COMM_WORLD, status,ierr)
        call mpi_recv(b(1,1), n*n, MPI_DOUBLE,0,2, MPI_COMM_WORLD, status,ierr)
    endif
!
    if(myrank.eq.0) then
        call mpi_send(a(1,1), n*n, MPI_DOUBLE,1,1,                   &
                    MPI_COMM_WORLD, ierr)
        call mpi_send(b(1,1), n*n, MPI_DOUBLE,1,2,                   &
                    MPI_COMM_WORLD, ierr)
    endif
    call mpi_barrier(MPI_COMM_WORLD,ierr)
!
    if (myrank == 0) then
        jstart = 1
        jend = n/2
    endif
!
    if (myrank == 1) then
        jstart = n/2+1
        jend = n
    endif
!
    do j=jstart, jend
        do k=1, n
            do i=1, n
                c(i,j) = c(i,j) + a(i,k)*b(k,j)
            end do
        end do
    end do
!
!collecting c elements
   if(myrank == 0) then
        call mpi_recv(c(1,n/2+1), n*n/2, MPI_DOUBLE,1,4,MPI_COMM_WORLD, status,ierr)
   endif
!
    if(myrank == 1) then
        call mpi_send(c(1,n/2+1), n*n/2, MPI_DOUBLE,0,4,MPI_COMM_WORLD, ierr)
    endif
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    time2 = MPI_Wtime()
    gflops = 2.0*n*n*n/(time2-time1)/1000**3
! check
    if(myrank == 0) then
        call random_number(d)
        i = int( d*n+1)
        call random_number(d)
        j = int( d*n+1)
        d = 0.d0
        do k=1, n
            d = d + a(i,k)*b(k,j)
        end do
    endif

    if(myrank == 0) then
        write(*,*) "[mm] Check on a random element:" , abs(d-c(i,j)), c(i,j)
        write(*,*) "[mm] Elapsed time ", time2-time1 ,' s'
        write(*,*) "[mm] Tot Gflops   ", gflops
        call MPI_SEND(gflops, 1, MPI_REAL8, 0, 1, MPI_COMM_MASTER, ierr)
    endif

    deallocate(a,b,c)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(myrank.eq.0) then
        call MPI_RECV(portname, 255, MPI_CHAR, 0, 2, MPI_COMM_MASTER, MPI_STATUS_IGNORE, ierr)
        write(*,'(a17,a100)') "[mm] portname = ", portname

        ! caution : use MPI_COMM_SELF, not MPI_COMM_WORLD for port connection
        ! create MPI_COMM_WORKER0 inter-communicator with pi.ex
        call MPI_COMM_CONNECT(portname, MPI_INFO_NULL, 0, MPI_COMM_SELF, MPI_COMM_WORKER0, ierr)
        call MPI_RECV(pi, 1, MPI_REAL8, 0, 1, MPI_COMM_WORKER0, MPI_STATUS_IGNORE, ierr)
        call MPI_SEND(gflops, 1, MPI_REAL8, 0, 1, MPI_COMM_WORKER0, ierr)
        write(*,*) "[mm] pi from pi.ex  ", pi

        call MPI_COMM_DISCONNECT(MPI_COMM_WORKER0, ierr)
    endif

    call MPI_FINALIZE(ierr)

end program mm
