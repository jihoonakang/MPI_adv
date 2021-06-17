program main

    use mpi_f08

    implicit none

    type(MPI_Comm)  :: comm_split, comm_split1, comm_split2
    type(MPI_Status):: status
    integer(kind=4) :: world_nprocs, world_myrank
    integer(kind=4) :: comm_split_nprocs, comm_split_myrank
    integer(kind=4) :: comm_split1_nprocs, comm_split1_myrank
    integer(kind=4) :: comm_split2_nprocs, comm_split2_myrank
    integer(kind=4) :: comm_split1_root_rank, comm_split2_root_rank

    integer(kind=4) :: color
    integer(kind=4) :: a, a_sum, b
    integer(kind=4) :: ierr, temp

    call MPI_Init()
    call MPI_Comm_size( MPI_COMM_WORLD, world_nprocs)
    call MPI_Comm_rank( MPI_COMM_WORLD, world_myrank)

    ! MPI_Reduce in global communicator
    a_sum = 0
    a = world_myrank
    call MPI_Reduce(a, a_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD)
    ! write(*,'(3(a,i5))') 'a = ', a, ', sum of a = ', a_sum, ', world_myrank = ',world_myrank

    ! Common split communicator 

    select case (mod(world_myrank,3))
    case (0)
        color = 0
    case default
        color = 1
    end select

    call MPI_Comm_split(MPI_COMM_WORLD, color, world_myrank, comm_split)
    call MPI_Comm_size( comm_split, comm_split_nprocs)
    call MPI_Comm_rank( comm_split, comm_split_myrank)
    ! write (*,'(4i5)') , world_nprocs, world_myrank, comm_split_nprocs, comm_split_myrank

    ! MPI_Reduce in common split communicator
    a_sum = 0
    a = world_myrank
    call MPI_Reduce(a, a_sum, 1, MPI_INTEGER, MPI_SUM, 0, comm_split)
    ! write(*,'(4(a,i5))') 'a = ', a, ', sum of a = ', a_sum, ', world_myrank = ',world_myrank,', comm_split_myrank = ',comm_split_myrank
 
    ! initialze nprocs and myrank in case of empty communicator
    comm_split1_myrank = MPI_PROC_NULL
    comm_split2_myrank = MPI_PROC_NULL
    comm_split1_nprocs = MPI_PROC_NULL
    comm_split2_nprocs = MPI_PROC_NULL

    ! individual split communicator
    select case (color)
    case (0)
        call MPI_Comm_dup(comm_split, comm_split1)
        call MPI_Comm_size( comm_split1, comm_split1_nprocs)
        call MPI_Comm_rank( comm_split1, comm_split1_myrank)
    case (1)
        call MPI_Comm_dup(comm_split, comm_split2)
        call MPI_Comm_size( comm_split2, comm_split2_nprocs)
        call MPI_Comm_rank( comm_split2, comm_split2_myrank)
    end select

    ! write (*,'(4i5)') , world_nprocs, world_myrank, comm_split2_nprocs, comm_split2_myrank
    ! write(*,'(4i5)'), world_nprocs, world_myrank, comm_split1_nprocs, comm_split1_myrank, comm_split2_nprocs, comm_split2_myrank
    a_sum = 0
    a = world_myrank

    select case (color)
    case (0)
        ! run program with MPI processes in comm_split1 communicator
        call MPI_Reduce(a, a_sum, 1, MPI_INTEGER, MPI_SUM, 0, comm_split1)
        write(*,'(4(a,i5))') 'a = ', a, ', sum of a = ', a_sum, ', world_myrank = ',world_myrank,', comm_split1_myrank = ',comm_split1_myrank
    case(1)
        ! run program with MPI processes in comm_split2 communicator
        call MPI_Reduce(a, a_sum, 1, MPI_INTEGER, MPI_SUM, 0, comm_split2)
        write(*,'(4(a,i5))') 'a = ', a, ', sum of a = ', a_sum, ', world_myrank = ',world_myrank,', comm_split2_myrank = ',comm_split2_myrank
    end select

    ! Finding world_myrank of root process of split communicator in MPI_COMM_WORLD
    temp = 0
    if(comm_split1_myrank.eq.0) then
        temp = world_myrank
    endif

    call MPI_Allreduce(temp, comm_split1_root_rank, 1 ,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
    write(*,*) comm_split1_root_rank, world_myrank

    temp = 0
    if(comm_split2_myrank.eq.0) then
        temp = world_myrank
    endif

    call MPI_Allreduce(temp, comm_split2_root_rank, 1 ,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
    write(*,*) comm_split2_root_rank, world_myrank

    ! p2p ommunication between root processes of split communicator in MPI_COMM_WORLD
    b = 0
    if(comm_split2_myrank.eq.0) then
        call MPI_Send(a_sum, 1, MPI_INTEGER, comm_split1_root_rank, 1, MPI_COMM_WORLD)
    endif

    if(comm_split1_myrank.eq.0) then
        call MPI_Recv(b, 1, MPI_INTEGER, comm_split2_root_rank, 1, MPI_COMM_WORLD, status)
    endif

    write(*,*) world_myrank, b

    call MPI_Comm_free(comm_split)
    select case (color)
    case (0)
        call MPI_Comm_free(comm_split1)
    case (1)
        call MPI_Comm_free(comm_split2)
    end select

    call MPI_Finalize()

end program