PROGRAM pi_integral
    include "mpif.h"

    integer, parameter:: num_steps=1000000000
    real(8) sum, step, x, pi, tsum;
    integer :: i, nprocs, myrank, ierr

    call MPI_INIT(ierr);
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

    sum=0.0
    step=1./dble(num_steps)

    call para_range(1, num_steps, nprocs, myrank, ista, iend)
    print*, "myrank =", myrank, ":", ista," ~ ", iend
    do i=ista, iend
        x = (i-0.5)*step
        sum = sum + 4.0/(1.0+x*x)
    enddo
    
    call MPI_REDUCE(sum, tsum, 1, MPI_REAL8, MPI_SUM, 0, & 
                    MPI_COMM_WORLD, ierr)

    if(myrank ==0) then
        pi = step*tsum
        print*, "numerical  pi = ", pi
        print*, "analytical pi = ", dacos(-1.d0)
        print*, " Error = ", dabs(dacos(-1.d0)-pi)
    endif

    call MPI_FINALIZE(ierr)
end 

SUBROUTINE para_range(n1, n2, nprocs, irank, ista, iend)
    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = MOD(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
    iend = ista + iwork1 - 1
    IF (iwork2 > irank) iend = iend + 1
END       