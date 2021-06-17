subroutine mpmd_comm(oldcomm, oldsize, oldrank, newcomm, newsize, newrank,color)

    use mpi
    implicit none

    integer,intent(in):: oldcomm, oldsize, oldrank
    integer,intent(inout):: newcomm, newsize, newrank, color
    integer :: ierr

    if(oldrank.LT.4) then
        color = 0
    else
        color = 1
    endif

    call MPI_COMM_SPLIT(oldcomm, color, oldrank, newcomm, ierr)
    call MPI_COMM_SIZE(newcomm,newsize,ierr)
    call MPI_COMM_RANK(newcomm,newrank,ierr)
    print *, '[common] myrank = ',oldrank, ' newrank = ',newrank

end subroutine mpmd_comm

subroutine mpmd_from_c0_to_c1_real8(sendrank, recvrank, myrank, a)

    use mpi 
    implicit none
    integer,intent(in):: sendrank, recvrank, myrank
    real(kind(1.d0)), intent(inout) :: a
    integer :: status(MPI_STATUS_SIZE), ierr

    if(myrank.eq.sendrank) then
        call MPI_Send(a, 1, MPI_REAL8, recvrank, 1, MPI_COMM_WORLD, ierr)
    endif


    if(myrank.eq.recvrank) then
        call MPI_Recv(a, 1, MPI_REAL8, sendrank, 1, MPI_COMM_WORLD, status, ierr)
    endif

end subroutine mpmd_from_c0_to_c1_real8
