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

inewtype = MPI.Datatype.Create_contiguous(MPI.INTEGER4, 2)
MPI.Datatype.Commit(inewtype)

comm.Bcast((ibuf, 3, inewtype), 0)

print(rank, ibuf)

MPI.Datatype.Free(inewtype)
