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

num_step = 1000000

dx = 1.0 / num_step

ista, iend = para_range(1, num_step, size, rank)

print('Rank = %d, (ista, iend) = (%d, %d)'%(rank, ista, iend))

sum = 0.0
for i in range(ista, iend+1) :
    x= (i - 0.5) * dx
    sum += 4.0/(1.0 + x*x)

tsum = comm.reduce(sum)

if rank == 0 :
    pi = dx * tsum
    print('Numerical pi = %f'%pi)
