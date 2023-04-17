################################################################################
# MPI python example
#
# Copyright 2023. Ji-Hoon Kang, All rights reserved.
# This project is released under the terms of the MIT License (see LICENSE )
################################################################################

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0 :
    ibuf = np.arange(1, 21, dtype = np.int32)
else :
    ibuf = np.zeros(20, dtype = np.int32)

inewtype = MPI.Datatype.Create_vector(MPI.INTEGER4, 1, 3, 5)
inewtype2 = MPI.Datatype.Create_resized(inewtype, 0, 5*4)
MPI.Datatype.Commit(inewtype)
MPI.Datatype.Commit(inewtype2)

comm.Bcast((ibuf, 4, inewtype2), 0)

print(rank, ibuf)

