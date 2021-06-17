! reference : https://www.hpc.cineca.it/content/mm-fortran-mpi

program matrix_matrix_prod
    use mpi
    implicit none
    integer :: n, size
    integer :: nprocs, ierr, myrank, status(MPI_STATUS_SIZE)
    real(kind(1.d0)), dimension(:,:), allocatable   :: a, b, c
    real(kind(1.d0)) :: d
    real(kind(1.d0)) :: time1, time2
    real(kind(1.d0)) :: pi

    integer            :: i, j, k, jstart, jend
    character(len=128) :: command
    character(len=80) :: arg

    integer             :: color
    integer             :: newcomm, newrank, newsize

!mpi stuff
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    call mpmd_comm(MPI_COMM_WORLD,nprocs,myrank,newcomm,newsize,newrank,color)

!reading matrix size
    if(newrank.eq.0) then
        write(*,*) "[mm] MPI tasks = ", newsize
        call get_command_argument(0,command)
        if (command_argument_count() /= 1) then
            write(0,*) '[mm] Usage:', trim(command), '   matrix size'
            call MPI_ABORT(newcomm, 911,ierr)
            call MPI_FINALIZE(ierr)
            stop
        else
            call get_command_argument(1,arg)
            read(arg,*) n
        endif
    endif
    write(*,*) "[mm] MPI ranks = ", newrank
    call MPI_BCAST(n,1,MPI_INTEGER,0,newcomm,ierr)
    call MPI_BARRIER(newcomm,ierr)
!a check
    if (n > 0 ) then
        if (newrank == 0 ) then
            write(*,*) '[mm] Matrix size is ', n
        endif
    else
        write(*,*) "[mm] Error, matrix size is ", n
        call MPI_FINALIZE(ierr)
        stop
    endif
!
    if(newsize /= 2) then
        if(newrank == 0) then
            write(*,*) "[mm] Error, newsize=", newsize, "is not 2!!!!"
            write(*,*) "[mm] This (stupid) code works only with 2 task!!!"
        endif
        call MPI_BARRIER(newcomm,ierr)
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

    if(newrank == 0) then
        call random_number(a)
        call random_number(b)
    endif
!
    time1 = MPI_Wtime()
!sending a and b elements
    if(newrank.eq.1) then
        call mpi_recv(a(1,1), n*n, MPI_DOUBLE,0,1, newcomm, status,ierr)
        call mpi_recv(b(1,1), n*n, MPI_DOUBLE,0,2, newcomm, status,ierr)
    endif
!
    if(newrank.eq.0) then
        call mpi_send(a(1,1), n*n, MPI_DOUBLE,1,1,                   &
                    newcomm, ierr)
        call mpi_send(b(1,1), n*n, MPI_DOUBLE,1,2,                   &
                    newcomm, ierr)
    endif
    call mpi_barrier(newcomm,ierr)
!
    if (newrank == 0) then
        jstart = 1
        jend = n/2
    endif
!
    if (newrank == 1) then
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
   if(newrank == 0) then
        call mpi_recv(c(1,n/2+1), n*n/2, MPI_DOUBLE,1,4,newcomm, status,ierr)
   endif
!
    if(newrank == 1) then
        call mpi_send(c(1,n/2+1), n*n/2, MPI_DOUBLE,0,4,newcomm, ierr)
    endif
    call mpi_barrier(newcomm,ierr)
    time2 = MPI_Wtime()
! check
    if(newrank == 0) then
        call random_number(d)
        i = int( d*n+1)
        call random_number(d)
        j = int( d*n+1)
        d = 0.d0
        do k=1, n
            d = d + a(i,k)*b(k,j)
        end do
    endif

    if(newrank == 0) then
        write(*,*) "[mm] Check on a random element:" , abs(d-c(i,j)), c(i,j)
        write(*,*) "[mm] Elapsed time ", time2-time1 ,' s'
        write(*,*) "[mm] Tot Gflops   ", 2.0*n*n*n/(time2-time1)/1000**3
    endif

    deallocate(a,b,c)

    call mpmd_from_c0_to_c1_real8(0,4,myrank,pi)
    if(newrank == 0) then
        write(*,*) "[mm] pi from pi.ex = ", pi
    endif
    if(newrank == 0) then
        write(*,*) "[mm] All done..."
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

end program matrix_matrix_prod
