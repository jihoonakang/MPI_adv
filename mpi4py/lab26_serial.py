import numpy as np
import matplotlib.pyplot as plt

global EPSILON
global MAX_ITER
global LX
global LY

def Jacobi_iter(N, M, psi_new, beta, beta_1) :

    EPSILON = 1e-8
    MAX_ITER = 1000

    psi_old = np.empty_like(psi_new)
    for iter in range(MAX_ITER) :
        error = 0.0
        psi_old[:, :] = psi_new[:, :]

        for i in range(1, N - 1) :
            for j in range(1, M - 1) :
                psi_new[i, j] = beta_1 * (psi_old[i, j + 1] + psi_old[i, j - 1] \
                              + beta * beta * (psi_old[i + 1, j] + psi_old[i - 1, j]))

        psi_new[:, M - 1] = psi_new[:, M - 2]

        for i in range(N) :
            for j in range(M) :
                error += (psi_new[i, j] - psi_old[i, j])**2
        error /= M * N
        error = error**(1/2) 

        if iter % 10 == 0 :
            print('Iteration = {0}, Error= {1}'.format(iter, error))
            # plt.clf()
            # plt.pcolormesh(psi_new, cmap=plt.cm.jet, vmin=0, vmax=100)
            # plt.colorbar()
            # plt.draw()
            # plt.pause(0.1)

        if error <= EPSILON :
            break


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

Jacobi_iter(N, M, psi_new, beta, beta_1)

#   for(i=0;i<divide;++i)
#     psi_new[N-1][i]=0.0; // bottom(left)
#   for(i=divide;i<M;++i)
#     psi_new[N-1][i]=100.0;  // bottom(right)
#   for(i=0;i<N;++i)
#     psi_new[i][0]=0.0;  // left wall
#   for(i=0;i<M;++i)
#     psi_new[0][i]=0.0;  // upper wall

# /* stream.c */
# #include <stdio.h>
# #include <math.h>
# const double LX=6.0, LY=4.0;  // length of domain along x-, y-direction
# const double EPSILON=1e-8;    // tolerance
# const long int MAX_ITER=1000000000;

# void Jacobi_iter(int N, int M, double psi_new[N][M],double beta,double beta_1);
# void pre_post(int N, int M, double psi_new[N][M], double dx, double dy);

# int main(void){
#   const int M=90*3+1, N=60*3+1; // # of grid points starting from 0 in x-, y-direction
#   double dx, dy, beta, beta_1;
#   int i, j, divide;
#   double psi_new[N][M];
#   double STime, ETime;

#   dx=LX/(M-1),  dy=LY/(N-1);
#   beta=(dx/dy);
#   beta_1=1.0/(2.0*(1.0+beta*beta));
#   for(i=0;i<N;++i)
#     for(j=0;j<M;++j)
#     {
#       psi_new[i][j]=0.0;
#     }
#  // boundary conditions
#   divide = (int)(M-1)*0.5;
#   for(i=0;i<divide;++i)
#     psi_new[N-1][i]=0.0; // bottom(left)
#   for(i=divide;i<M;++i)
#     psi_new[N-1][i]=100.0;  // bottom(right)
#   for(i=0;i<N;++i)
#     psi_new[i][0]=0.0;  // left wall
#   for(i=0;i<M;++i)
#     psi_new[0][i]=0.0;  // upper wall
#   // Jacobi_ieration
#   Jacobi_iter(N,M,psi_new,beta,beta_1);

#   // write mesh & post-processihng file
#   pre_post(N, M, psi_new, dx, dy);

#   return 0;
# }
# /**********************************************************************************/

# void Jacobi_iter(int N, int M,double psi_new[N][M],double beta,double beta_1)
# {
#   long int iter;
#   double error=1.0;
#   int i, j;
#   double psi_old[N][M];


#   for(iter=1;iter<MAX_ITER;++iter)
#   {
#     error=0.0;
#     for(i=0;i<N;++i)
#       for(j=0;j<M;++j)
#         psi_old[i][j]=psi_new[i][j];

#     for(i=1;i<N-1;++i)
#       for(j=1;j<M-1;++j)
#         psi_new[i][j]=beta_1*(psi_old[i][j+1]+psi_old[i][j-1]+
#                               beta*beta*(psi_old[i+1][j]+psi_old[i-1][j]));

#  // Right Neumann Boundary Condition
#     for(i=0;i<N;++i)
#       psi_new[i][M-1]=psi_new[i][M-2];

#     for(i=0;i<N;++i)
#       for(j=0;j<M;++j)
#         error += (psi_new[i][j]-psi_old[i][j])*(psi_new[i][j]-psi_old[i][j]);
#     error = sqrt(error/(M*N));
#     if(iter%1000==0) printf("Iteration = %ld, Error= %lg\n",iter,error);

#     if(error<=EPSILON) break;
#   }
#   printf("Iteration = %ld, Error= %lg\n",iter,error);
# }


# /* stream.c */
# #include <stdio.h>
# #include <math.h>
# const double LX=6.0, LY=4.0;  // length of domain along x-, y-direction
# const double EPSILON=1e-8;    // tolerance
# const long int MAX_ITER=1000000000;

# void Jacobi_iter(int N, int M, double psi_new[N][M],double beta,double beta_1);
# void pre_post(int N, int M, double psi_new[N][M], double dx, double dy);

# int main(void){
#   const int M=90*3+1, N=60*3+1; // # of grid points starting from 0 in x-, y-direction
#   double dx, dy, beta, beta_1;
#   int i, j, divide;
#   double psi_new[N][M];
#   double STime, ETime;

#   dx=LX/(M-1),  dy=LY/(N-1);
#   beta=(dx/dy);
#   beta_1=1.0/(2.0*(1.0+beta*beta));
#   for(i=0;i<N;++i)
#     for(j=0;j<M;++j)
#     {
#       psi_new[i][j]=0.0;
#     }
#  // boundary conditions
#   divide = (int)(M-1)*0.5;
#   for(i=0;i<divide;++i)
#     psi_new[N-1][i]=0.0; // bottom(left)
#   for(i=divide;i<M;++i)
#     psi_new[N-1][i]=100.0;  // bottom(right)
#   for(i=0;i<N;++i)
#     psi_new[i][0]=0.0;  // left wall
#   for(i=0;i<M;++i)
#     psi_new[0][i]=0.0;  // upper wall
#   // Jacobi_ieration
#   Jacobi_iter(N,M,psi_new,beta,beta_1);

#   // write mesh & post-processihng file
#   pre_post(N, M, psi_new, dx, dy);

#   return 0;
# }
# /**********************************************************************************/

# /***********************************************************************/
# void pre_post(int N, int M, double  psi_new[N][M], double dx, double dy)
# {
# 	0. https://gmsh.info/
# 	1. Gmsh mesh file 생성
# 	2. Gmsh post-proces file 생성
#   …
# }
