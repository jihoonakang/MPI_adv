! reference : https://www.hpc.cineca.it/content/pi-fortran-mpi
!---------------------------------------------c
!  Exercise: Pi                               c
!                                             c
!  Compute the value of PI using the integral c
!  pi = 4* int 1/(1+x*x)    x in [0-1]        c
!                                             c
!  The integral is approximated by a sum of   c
!  n interval.                                c
!                                             c
!  The approximation to the integral in each  c
!  interval is: (1/n)*4/(1+x*x).              c
!---------------------------------------------c
program pigreco
    use mpi
    implicit none

    integer(selected_int_kind(18)) :: i, istart, iend
    integer(selected_int_kind(18)), parameter :: intervals=1e7
    integer:: ierr, nprocs, myrank
    real(kind(1.d0)) :: dx,sum,x, total_sum
    real(kind(1.d0)) :: f,pi
    real(kind(1.d0)), parameter :: PI25DT = acos(-1.d0)
    real(kind(1.d0)) :: time1, time2

    integer             :: color
    integer             :: newcomm, newrank, newsize
    integer             :: MPI_COMM_MASTER

!mpi stuff
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    call MPI_Comm_get_parent(MPI_COMM_MASTER,ierr)

    call mpmd_comm(MPI_COMM_WORLD,nprocs,myrank,newcomm,newsize,newrank,color)

    if(newrank == 0) then
       write(*,*) "[pi] MPI tasks = ", newsize
       write(*,*) "[pi] Number of intervals     = ", intervals
    endif
	write(*,*) "[pi] MPI ranks = ", myrank
    if(mod(intervals,newsize) /= 0) then
       if(newrank == 0) then
          write(*,*) "[pi] The number of process must divide", intervals, "exactly."
       endif
       call MPI_BARRIER(newcomm,ierr)
       call MPI_FINALIZE(ierr)
       stop
    endif

    sum=0.d0
    dx=1.d0/intervals
    time1 = MPI_WTIME()
    istart = (intervals/newsize)*newrank + 1
    iend   = (intervals/newsize)*(newrank+1)
    sum = 0.0;
    total_sum = 0.0;
    do i=iend, istart, -1
        x=dx*(i-0.5d0)
        f=4.d0/(1.d0+x*x)
        sum=sum+f
    end do
    call MPI_Reduce(sum,total_sum,1,MPI_DOUBLE,MPI_SUM,0,newcomm,ierr)
    time2 = MPI_WTIME()
    pi=dx*total_sum
    call MPI_BARRIER(newcomm,ierr)
    if(newrank == 0) then
       PRINT '(a19,2x,f30.25)','[pi] Computed PI =', pi
       PRINT '(a19,2x,f30.25)','[pi] The True PI =', PI25DT
       PRINT '(a19,2x,f30.25)','[pi] Error        ', PI25DT-pi
       PRINT *, '[pi] Elapsed time ', time2-time1 ,' s'
       call MPI_SEND(pi, 1, MPI_REAL8, 0, 1, MPI_COMM_MASTER, ierr)
    endif

    call mpmd_from_c0_to_c1_integer(0,4,myrank,pi)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

contains

end program
