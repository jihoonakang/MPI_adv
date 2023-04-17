################################################################################
# MPI python example
#
# Copyright 2023. Ji-Hoon Kang, All rights reserved.
# This project is released under the terms of the MIT License (see LICENSE )
################################################################################

def para_range(n1, n2, size, rank) :
    iwork = divmod((n2 - n1 + 1), size)
    ista = rank * iwork[0] + n1 + min(rank, iwork[1])
    iend = ista + iwork[0] - 1
    if iwork[1] > rank :
        iend = iend + 1

    return ista, iend

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()

n = 11

a = np.zeros(n, dtype = np.int32)
b = np.zeros(n, dtype = np.int32)

ista, iend = para_range(0, n - 1, size, rank)

for i in range(ista, iend+1) :
    b[i] = i + 1

ista1 = ista; iend1 = iend

if rank == 0 :
    ista1 = 1
if rank == size - 1 :
    iend1 = n - 2

inext = rank + 1; iprev = rank - 1

if rank == size - 1 :
    inext = MPI.PROC_NULL
if rank == 0 :
    iprev = MPI.PROC_NULL

req_i1 = comm.Isend(b[iend:iend+1], inext, 11)
req_i2 = comm.Isend(b[ista:ista+1], iprev, 12)
req_r1 = comm.Irecv(b[ista-1: ista], iprev, 11)
req_r2 = comm.Irecv(b[iend+1: iend+2], inext, 12)

MPI.Request.Wait(req_i1)
MPI.Request.Wait(req_i2)
MPI.Request.Wait(req_r1)
MPI.Request.Wait(req_r2)

for i in range(ista1, iend1+1) :
    a[i] = b[i-1] + b[i+1]

for i in range(size) :
    if i == rank :
        print(rank)
        print(b)
        print(a)
    comm.Barrier()

