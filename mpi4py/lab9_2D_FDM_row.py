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

m = 6; n = 9

a = np.zeros((m, n), dtype = np.float64)
b = np.zeros((m, n), dtype = np.float64)

ista, iend = para_range(0, m - 1, size, rank)

for i in range(ista, iend+1) :
    for j in range(n) :
        a[i,j] = i + 1 + 10 * j

ista1 = ista; iend1 = iend

if rank == 0 :
    ista1 = 1
if rank == size - 1 :
    iend1 = m - 2

inext = rank + 1; iprev = rank - 1

if rank == size - 1 :
    inext = MPI.PROC_NULL
if rank == 0 :
    iprev = MPI.PROC_NULL

if rank != size-1 :
    works1 = a[:,iend].copy()
else :
    works1 = np.zeros(m, dtype = np.float64)

if(rank != 0) : 
    works2 = a[:,iend].copy()
else :
    works2= np.zeros(m, dtype = np.float64)

workr1 = np.zeros(m, dtype = np.float64)
workr2 = np.zeros(m, dtype = np.float64)

if(rank != 0) : 
    req_i2 = comm.Isend(a[ista,:], iprev, 12)
    req_r1 = comm.Irecv(a[ista-1,:], iprev, 11)
if rank != size-1 :
    req_i1 = comm.Isend(a[iend,:], inext, 11)
    req_r2 = comm.Irecv(a[iend+1,:], inext, 12)

if(rank != 0) : 
    MPI.Request.Wait(req_i2)
    MPI.Request.Wait(req_r1)
if rank != size-1 :
    MPI.Request.Wait(req_i1)
    MPI.Request.Wait(req_r2)

for i in range(ista1, iend1+1):
    for j in range(1, n-1) :
        b[i][j] = a[i-1][j] + a[i][j-1] + a[i][j+1] + a[i+1][j]

for i in range(size) :
    if i == rank :
        print(rank)
        print(a)
        print(b)
    comm.Barrier()
