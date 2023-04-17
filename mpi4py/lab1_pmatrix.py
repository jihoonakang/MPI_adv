################################################################################
# MPI python example
#
# Copyright 2023. Ji-Hoon Kang, All rights reserved.
# This project is released under the terms of the MIT License (see LICENSE )
################################################################################

from mpi4py import MPI
import numpy as np
import random as rd

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()

NP = 5

matrixA = np.zeros((NP, NP), dtype = np.int32)
matrixB = np.zeros((NP, NP), dtype = np.int32)
matrixC = np.zeros((NP, NP), dtype = np.int32)
matrixT = np.zeros((NP, NP), dtype = np.int32)

if rank == 0 :
    for i in range(NP) :
        for j in range(NP) :
                matrixA[i][j] = rd.randrange(1, 10)
                matrixB[i][j] = rd.randrange(1, 10)
    print('Rank = %d'%rank)
    print(matrixA)

row_new_type = MPI.Datatype.Create_contiguous(MPI.INTEGER4, NP)
MPI.Datatype.Commit(row_new_type)

comm.Bcast((matrixA, 5, row_new_type), 0)
comm.Bcast((matrixB, 5, row_new_type), 0)

MPI.Datatype.Free(row_new_type)

for i in range(NP) :
    for j in range(NP) :
        matrixC[rank][i] = matrixC[rank][i] + matrixA[rank][j] * matrixB[j][i]

comm.Reduce(matrixC, matrixT, MPI.SUM, 0)

if rank == 0 :
    print('Rank = %d, parallel solution ='%rank)
    print(matrixT)
    print('Rank = %d, sequential solution ='%rank)
    print(matrixA@matrixB)


