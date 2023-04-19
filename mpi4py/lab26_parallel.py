from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

global EPSILON
global MAX_ITER
global LX
global LY

def para_range(n1, n2, size, rank) :
    iwork = divmod((n2 - n1 + 1), size)
    ista = rank * iwork[0] + n1 + min(rank, iwork[1])
    iend = ista + iwork[0] - 1
    if iwork[1] > rank :
        iend = iend + 1
    return ista, iend

def Jacobi_iter(N, M, psi_new, beta, beta_1, rank, size) :

    EPSILON = 1e-8
    MAX_ITER = 1000

    ista, iend = para_range(0, N - 1, size, rank)

    ista_p = ista
    iend_m = iend
    iprev = rank - 1
    inext = rank + 1

    if rank == 0 :
        ista_p = ista + 1
        iprev = MPI.PROC_NULL

    if rank == size - 1 :
        iend_m = iend - 1
        inext = MPI.PROC_NULL

    psi_old = np.empty_like(psi_new)

    for iter in range(MAX_ITER) :
        error_loc = 0.0
        psi_old[ista:iend+1, :] = psi_new[ista:iend+1, :]

        req_list = []
        if inext != MPI.PROC_NULL :
            req_send1 = comm.Isend((psi_old[iend,:], MPI.DOUBLE), inext, 1)
            req_recv2 = comm.Irecv((psi_old[iend+1,:], MPI.DOUBLE), inext, 2)
            req_list.append(req_send1)
            req_list.append(req_recv2)
        if iprev != MPI.PROC_NULL :
            req_recv1 = comm.Irecv((psi_old[ista-1,:], MPI.DOUBLE), iprev, 1)
            req_send2 = comm.Isend((psi_old[ista,:], MPI.DOUBLE), iprev, 2)
            req_list.append(req_send2)
            req_list.append(req_recv1)

        MPI.Request.Waitall(req_list)

        for i in range(ista_p, iend_m + 1) :
            for j in range(1, M - 1) :
                psi_new[i, j] = beta_1 * (psi_old[i, j + 1] + psi_old[i, j - 1] \
                              + beta * beta * (psi_old[i + 1, j] + psi_old[i - 1, j]))

        psi_new[ista:iend+1, M - 1] = psi_new[ista:iend+1, M - 2]

        for i in range(ista, iend + 1) :
            for j in range(M) :
                error_loc += (psi_new[i, j] - psi_old[i, j])**2
        error = comm.allreduce(error_loc, MPI.SUM)
        error /= M * N
        error = error**(1/2)


        if iter % 10 == 0 :

            ircnt = comm.gather((iend - ista + 1) * M)

            if rank == 0 :
                comm.Gatherv(MPI.IN_PLACE, (psi_new, ircnt), 0)
            else :
                comm.Gatherv(psi_new[ista:iend+1, :], 0)

            if rank == 0 :
                print('Iteration = {0}, Error= {1}'.format(iter, error))
                plt.clf()
                plt.pcolormesh(psi_new, cmap=plt.cm.jet, vmin=0, vmax=100)
                plt.colorbar()
                plt.draw()
                plt.pause(0.1)

        if error <= EPSILON :
            break

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

LX = 6.0
LY = 4.0

M = 90 * 1 + 1
N = 60 * 1 + 1

psi_new = np.zeros((N, M), dtype = np.double)

dx = LX / (M - 1)
dy = LY / (N - 1)
beta = dx / dy
beta_1 = 1.0 / (2.0 * (1.0 + beta * beta))

divide = int((M - 1) / 2)

for i in range(divide) :
    psi_new[N - 1, i] = 0.0
for i in range(divide, M) :
    psi_new[N - 1, i] = 100.0
for i in range(M) :
    psi_new[0, i] = 0.0

for i in range(N) :
    psi_new[i, 0] = 0.0

Jacobi_iter(N, M, psi_new, beta, beta_1, rank, size)
