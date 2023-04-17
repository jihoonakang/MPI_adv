################################################################################
# MPI python example
#
# Copyright 2023. Ji-Hoon Kang, All rights reserved.
# This project is released under the terms of the MIT License (see LICENSE )
################################################################################

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD

rank = comm.Get_rank()
size = comm.Get_size()

bufsize = 100
buf = np.zeros(bufsize, dtype = np.int32)

amode = MPI.MODE_RDONLY
f = MPI.File.Open(comm = MPI.COMM_SELF, filename = 'testfile', amode = amode)
disp = rank * bufsize * MPI.INTEGER4.Get_size()
MPI.File.Set_view(f, disp, MPI.INTEGER4, MPI.INTEGER4)
MPI.File.Read(f, buf)
MPI.File.Close(f)
print('Rank %d load data from testfile \n'%rank, buf)
