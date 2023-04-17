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

SCOPE = 1000000

mycount = 0
for i in range(SCOPE) :
    x = np.random.rand()
    y = np.random.rand()
    z = (x*x + y*y)**(0.5)
    if z < 1 :
        mycount += 1

count = comm.reduce(mycount)

if rank == 0 :
    print('Rank : %d, Count = %d, Pi = %f'%(rank,count,count/SCOPE/size*4))
