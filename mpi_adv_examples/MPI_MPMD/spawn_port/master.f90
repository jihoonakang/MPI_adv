program master

    use mpi 
    implicit none


    integer :: info(2), worker_np(2)
    integer :: intercomm, ierr, errcode
    character*25 :: worker(2)
    character*25 :: array_of_argv(2,1), argv(2)
    real(kind(1.d0)) :: gflops, pi

    integer :: MPI_COMM_WORKER(2)
    character(len=255) :: portname

    call MPI_INIT(ierr)
    call MPI_Info_create(info(1),ierr)
    call MPI_Info_set(info(1),"add-host","jihoon",ierr)
    call MPI_Info_create(info(2),ierr)
    call MPI_Info_set(info(2),"add-host","jihoon",ierr)

    ! info(1) = MPI_INFO_NULL
    ! info(2) = MPI_INFO_NULL

    worker(1) = './pi.ex'
    worker(2) = './mm.ex'

    worker_np(1) = 4
    worker_np(2) = 2

    argv(1) = ''
    argv(2) = '' ! Some MPI implementation require this

    call MPI_Comm_spawn(worker(1), argv, worker_np(1), info(1), 0, MPI_COMM_WORLD, MPI_COMM_WORKER(1), MPI_ERRCODES_IGNORE, ierr)
    call MPI_RECV(pi, 1, MPI_REAL8, 0, 1, MPI_COMM_WORKER(1), MPI_STATUS_IGNORE, ierr)
    call MPI_RECV(portname, 255, MPI_CHAR, 0, 2, MPI_COMM_WORKER(1), MPI_STATUS_IGNORE, ierr)

    write(*,*) "[master] pi   = ", pi
    write(*,'(a21,a100)') "[master] portname = ", portname

    argv(1) = '20'
    argv(2) = '' ! Some MPI implementation require this
    call MPI_Comm_spawn(worker(2), argv, worker_np(2), info(2), 0, MPI_COMM_WORLD, MPI_COMM_WORKER(2), MPI_ERRCODES_IGNORE, ierr)
    call MPI_SEND(portname, 255, MPI_CHAR, 0, 2, MPI_COMM_WORKER(2), ierr)
    call MPI_RECV(gflops, 1, MPI_REAL8, 0, 1, MPI_COMM_WORKER(2), MPI_STATUS_IGNORE, ierr)
    write(*,*) "[master] Tot Gflops   ", gflops

    call MPI_Finalize(ierr)

end program master

