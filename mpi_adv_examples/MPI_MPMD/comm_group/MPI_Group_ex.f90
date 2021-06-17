PROGRAM group

    USE MPI_F08
    IMPLICIT NONE

    INTEGER(KIND=4) :: world_rank, world_size, n, ranks(6)
    INTEGER(KIND=4) :: prime_rank = MPI_UNDEFINED, prime_size = MPI_UNDEFINED
    TYPE(MPI_GROUP) :: world_group, prime_group
    TYPE(MPI_COMM)  :: prime_comm
    INTEGER(KIND=4) :: prime_me

    CALL MPI_INIT() 

    ! Get the rank and size in the original communicator
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size) 
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, world_rank) 

    ! Get the group of processes in MPI_COMM_WORLD
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_group)

    n = 6
    ranks = (/2, 3, 5, 7, 11, 13/)

    ! Construct a group containing all of the prime ranks in world_group
    CALL MPI_GROUP_INCL(world_group, n, ranks, prime_group)
    CALL MPI_GROUP_RANK(prime_group, prime_me)
    ! If this rank isn't in the new group, me will be MPI_UNDEFINED.

    ! Create a new communicator based on the group
    ! CALL MPI_COMM_CREATE(MPI_COMM_WORLD, prime_group, prime_comm)
    IF (prime_me.NE.MPI_UNDEFINED) THEN
        CALL MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD, prime_group, 1, prime_comm)
    ENDIF

    ! If this rank isn't in the new communicator, it will be
    ! MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
    ! MPI_Comm_size is erroneous

    IF (prime_me.NE.MPI_UNDEFINED) THEN
        CALL MPI_COMM_RANK(prime_comm, prime_rank)
        CALL MPI_COMM_SIZE(prime_comm, prime_size)
    ENDIF
    PRINT '(A,2I7,A,2I7)', 'WORLD RANK/SIZE: ', world_rank, world_size, ' PRIME RANK/SIZE: ', prime_rank, prime_size
    IF (prime_me.NE.MPI_UNDEFINED) THEN
        CALL MPI_COMM_FREE(prime_comm)
    ENDIF
    IF (prime_me.NE.MPI_UNDEFINED) THEN
        CALL MPI_Group_free(prime_group)
    ENDIF
    
    CALL MPI_Group_free(world_group)
    CALL MPI_Finalize()
END

