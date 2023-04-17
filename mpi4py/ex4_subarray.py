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

ibuf = np.zeros((6, 7), dtype = int)

if rank == 0 :
    ibuf[0:6, 0:7] = np.arange(1, 8, dtype = 'd')

inewtype = MPI.Datatype.Create_subarray(MPI.INTEGER8, (6, 7), (2, 5), (1, 1), MPI.ORDER_C)
MPI.Datatype.Commit(inewtype)

comm.Bcast((ibuf, inewtype), 0)

print(rank, ibuf)

