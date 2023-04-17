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

jsta, jend = para_range(0, n - 1, size, rank)

for i in range(m) :
    for j in range(jsta, jend+1) :
        a[i,j] = i + 1 + 10 * j

jsta1 = jsta; jend1 = jend

if rank == 0 :
    jsta1 = 1
if rank == size - 1 :
    jend1 = n - 2

inext = rank + 1; iprev = rank - 1

if rank == size - 1 :
    inext = MPI.PROC_NULL
if rank == 0 :
    iprev = MPI.PROC_NULL

if rank != size-1 :
    works1 = a[:,jend].copy()
else :
    works1 = np.zeros(m, dtype = np.float64)

if(rank != 0) : 
    works2 = a[:,jsta].copy()
else :
    works2= np.zeros(m, dtype = np.float64)

workr1 = np.zeros(m, dtype = np.float64)
workr2 = np.zeros(m, dtype = np.float64)

req_i1 = comm.Isend(works1, inext, 11)
req_i2 = comm.Isend(works2, iprev, 12)
req_r1 = comm.Irecv(workr1, iprev, 11)
req_r2 = comm.Irecv(workr2, inext, 12)

MPI.Request.Wait(req_i1)
MPI.Request.Wait(req_i2)
MPI.Request.Wait(req_r1)
MPI.Request.Wait(req_r2)

if rank != 0 :
    a[:,jsta-1] = workr1[:]
if rank != size - 1 : 
    a[:,jend+1] = workr2[:]

for i in range(1, m-1):
    for j in range(jsta1, jend1+1) :
        b[i][j] = a[i-1][j] + a[i][j-1] + a[i][j+1] + a[i+1][j]

for i in range(size) :
    if i == rank :
        print(rank)
        print(a)
        print(b)
    comm.Barrier()
