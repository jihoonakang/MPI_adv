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

N = 20
a = np.zeros((N, N), dtype = np.float64)
disp = np.zeros(N, dtype = np.int32)
blocklen = np.zeros(N, dtype = np.int32)

if rank == 0 :
    for i in range(N) :
        for j in range(N) :
            a[i][j]= np.double(i + j + 1)

for i in range(N) :
    disp[i] = N * i + i
    blocklen[i] = N - i

upper = MPI.Datatype.Create_indexed(MPI.DOUBLE, blocklen, disp)
MPI.Datatype.Commit(upper)

if rank == 0 :
    comm.Send((a, 1, upper), 1, 99)
else :
    comm.Recv((a, 1, upper), 0, 99)

MPI.Datatype.Free(upper)

if rank == 1 :
    print('myrank = %d, buf = '%rank)
    for i in range(N) :
        for j in range(N) :
            print('%5.2f'%a[i][j], end=' ')
        print()
    
