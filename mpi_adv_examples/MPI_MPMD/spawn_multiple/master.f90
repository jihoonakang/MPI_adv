program master

    use mpi 
    implicit none


    integer :: info(2), worker_np(2)
    integer :: intercomm, ierr, errcode
    character*25 :: worker(2)
    character*25 :: array_of_argv(2,1), argv(1)
    real(kind(1.d0)) :: gflops, pi

    integer :: MPI_COMM_WORKER

    call MPI_INIT(ierr)
    call MPI_Info_create(info(1),ierr)
    call MPI_Info_set(info(1),"add-host","jihoon",ierr)
    call MPI_Info_create(info(2),ierr)
    call MPI_Info_set(info(2),"add-host","jihoon",ierr)

    worker(1) = './pi.ex'
    worker(2) = './mm.ex'

    array_of_argv(1,1) = ''
    array_of_argv(2,1) = '20'

    worker_np(1) = 4
    worker_np(2) = 2

    call MPI_Comm_spawn_multiple(2,worker,array_of_argv, worker_np, info, 0, MPI_COMM_WORLD,MPI_COMM_WORKER,MPI_ERRCODES_IGNORE,ierr)
    call MPI_RECV(pi, 1, MPI_REAL8, 0, 1, MPI_COMM_WORKER, MPI_STATUS_IGNORE, ierr)
    write(*,*) "[master] pi   ", pi

    call MPI_RECV(gflops, 1, MPI_REAL8, 4, 1, MPI_COMM_WORKER, MPI_STATUS_IGNORE, ierr)
    write(*,*) "[master] Tot Gflops   ", gflops

    call MPI_Finalize(ierr)

end program master

