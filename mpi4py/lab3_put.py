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

a = np.array(pow(10, rank), dtype = np.int32)

print('Rank = %d, a = %d'%(rank, a))

win = MPI.Win.Create(a, comm = comm)

win.Fence()

if rank == 0 :
    win.Put(a, target_rank = 1)
    win.Put(a, target_rank = 2)

win.Fence()

print('Rank after put = %d, a = %d'%(rank,a))
