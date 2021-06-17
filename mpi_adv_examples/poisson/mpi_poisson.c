#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

typedef struct mympi {
    int nprocs;
    int myrank;
    int nx_mpi;
    int ny_mpi;
    int mpisize_x;
    int mpisize_y;
    int mpirank_x;
    int mpirank_y;
    int w_rank;
    int e_rank;
    int n_rank;
    int s_rank;
    MPI_Datatype vector_type;
    MPI_Comm comm_x;
    MPI_Comm comm_y;

} MYMPI;

void mpi_setup(int nx, int ny, MYMPI *mpi_info);
void send_east(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_west(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_north(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_south(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);

void RB_gauss_seidel(double dx, double dy, int nx, int ny, double *u_solve, double *rhs, int maxiteration, double tolerance, MYMPI *mpi_info)
{
    double *u_new;
    double error_sum, u_sum;
    double error_sum_global, u_sum_global;
    double ae,aw,as,an,ap;
    double dxsqi,dysqi;
    double uij_old,uij_new;
    int i,j,iter,js,walk;

    dxsqi=1.0/dx/dx;
    dysqi=1.0/dy/dy;

    for(iter=0;iter<maxiteration;iter++)
    {
        error_sum = 0.0;
        error_sum_global = 0.0;
        u_sum = 0.0;
        u_sum_global = 0.0;

        for(walk=1;walk<3;walk++)
        {
            if(mpi_info->w_rank < 0 )
            {
                for(j=1;j<=ny;j++)
                {
                    u_solve[0*(ny+2)+j]=-u_solve[1*(ny+2)+j];
                }
            }
            if(mpi_info->e_rank < 0 )
            {
                for(j=1;j<=ny;j++)
                {
                    u_solve[(nx+1)*(ny+2)+j]=-u_solve[nx*(ny+2)+j];
                }
            }
            send_east(u_solve, nx, ny, mpi_info);
            send_west(u_solve, nx, ny, mpi_info);
            send_north(u_solve, nx, ny, mpi_info);
            send_south(u_solve, nx, ny, mpi_info);

            js = 3 - walk;
            for(i=1;i<=nx;i++)
            {
                js=3-js;
                for(j=js;j<=ny;j+=2)
                {
                    aw=dxsqi*u_solve[(i-1)*(ny+2)+j];
                    ae=dxsqi*u_solve[(i+1)*(ny+2)+j];
                    as=dysqi*u_solve[i*(ny+2)+j-1];
                    an=dysqi*u_solve[i*(ny+2)+j+1];

                    ap = 2.0*(dxsqi+dysqi);

                    uij_old = u_solve[i*(ny+2)+j];
                    uij_new = (aw+ae+as+an-rhs[i*(ny+2)+j]) /ap;

                    error_sum += fabs(uij_new-uij_old);
                    u_sum += fabs(uij_new);

                    u_solve[i*(ny+2)+j] = uij_new;
                }
            }
        }

        MPI_Allreduce(&error_sum,&error_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&u_sum,&u_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(mpi_info->myrank ==0 && iter%100==0)
        {
            printf("%5d th iteration: error=%e\n",iter,error_sum_global/u_sum_global);
        }
        if(error_sum_global/u_sum_global<tolerance) break;
    }
    if(mpi_info->myrank ==0)
    printf("Iteration ends in %d th step\n",iter);
}

void mpi_setup(int nx, int ny, MYMPI *mpi_info)
{
    int dimsize[2],periods[2],coords[2],remain[2],reorder;
    MPI_Comm comm_cart;

    mpi_info->mpisize_x=2;
    mpi_info->mpisize_y=4;
    mpi_info->nx_mpi=nx/mpi_info->mpisize_x;
    mpi_info->ny_mpi=ny/mpi_info->mpisize_y;

    dimsize[0]=mpi_info->mpisize_x;
    dimsize[1]=mpi_info->mpisize_y;
    periods[0]=0;
    periods[1]=1;
    reorder=0;
    MPI_Cart_create(MPI_COMM_WORLD,2,dimsize,periods,reorder,&comm_cart);
    MPI_Cart_coords(comm_cart,mpi_info->myrank,2,coords);
    mpi_info->mpirank_x=coords[0];
    mpi_info->mpirank_y=coords[1];

    MPI_Cart_shift(comm_cart,0,1,&mpi_info->w_rank,&mpi_info->e_rank);
    MPI_Cart_shift(comm_cart,1,1,&mpi_info->s_rank,&mpi_info->n_rank);

    MPI_Type_vector(mpi_info->nx_mpi,1,mpi_info->ny_mpi+2,MPI_DOUBLE,&mpi_info->vector_type);
    MPI_Type_commit(&mpi_info->vector_type);

    remain[0]=0;
    remain[1]=1;
    MPI_Cart_sub(comm_cart,remain,&mpi_info->comm_y);
    MPI_Cart_shift(mpi_info->comm_y,0,1,&mpi_info->s_rank,&mpi_info->n_rank);

    remain[0]=1;
    remain[1]=0;
    MPI_Cart_sub(comm_cart,remain,&mpi_info->comm_x);
    MPI_Cart_shift(mpi_info->comm_x,0,1,&mpi_info->w_rank,&mpi_info->e_rank);
}
void send_east(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int j;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    if(mpi_info->e_rank >= 0)
    {
        MPI_Isend(&u[nx_mpi*(ny_mpi+2)+1],ny_mpi,MPI_DOUBLE,mpi_info->e_rank,101,mpi_info->comm_x,&req1);
    }
    if(mpi_info->w_rank >= 0)
    {
        MPI_Irecv(&u[0*(ny_mpi+2)+1],ny_mpi,MPI_DOUBLE,mpi_info->w_rank,101,mpi_info->comm_x,&req2);
    }
    if(mpi_info->e_rank >= 0)
    {
        MPI_Wait(&req1,&status1);
    }
    if(mpi_info->w_rank >= 0)
    {
        MPI_Wait(&req2,&status2);
    }
}

void send_west(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int j;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    if(mpi_info->w_rank >= 0)
    {
        MPI_Isend(&u[1*(ny_mpi+2)+1],ny_mpi,MPI_DOUBLE,mpi_info->w_rank,102,mpi_info->comm_x,&req1);
    }
    if(mpi_info->e_rank >= 0)
    {
        MPI_Irecv(&u[(nx_mpi+1)*(ny_mpi+2)+1],ny_mpi,MPI_DOUBLE,mpi_info->e_rank,102,mpi_info->comm_x,&req2);
    }
    if(mpi_info->w_rank >= 0)
    {
        MPI_Wait(&req1,&status1);
    }
    if(mpi_info->e_rank >= 0)
    {
        MPI_Wait(&req2,&status2);
    }
}

void send_north(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    MPI_Isend(&u[(ny_mpi+2)+ny_mpi],1,mpi_info->vector_type,mpi_info->n_rank,103,mpi_info->comm_y,&req1);
    MPI_Irecv(&u[(ny_mpi+2)],1,mpi_info->vector_type,mpi_info->s_rank,103,mpi_info->comm_y,&req2);

    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);
}

void send_south(double *u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    MPI_Isend(&u[(ny_mpi+2)+1],1,mpi_info->vector_type,mpi_info->s_rank,104,mpi_info->comm_y,&req1);
    MPI_Irecv(&u[(ny_mpi+2)+ny_mpi+1],1,mpi_info->vector_type,mpi_info->n_rank,104,mpi_info->comm_y,&req2);
    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);
}

int main(int argc, char **argv)
{
    int nx = 256, ny = 256;
    int maxiter = 100000;
    double tol = 1.0e-8;
    double length_x = 1.0, length_y = 1.0;
    double PI = atan(1.0)*4.0;

    int grid_size,grid_size_mpi;
    double dx, dy, x_val, y_val;
    double *pos_x, *pos_y; double *u_exact, *u_solve, *rhs; 
    int i,j;

    MYMPI mpi_info;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_info.nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_info.myrank);

    mpi_setup(nx,ny,&mpi_info);

    printf("%d %d %d %d %d %d %d %d %d %d\n",mpi_info.nprocs,mpi_info.myrank,mpi_info.mpisize_x,mpi_info.mpisize_y,mpi_info.mpirank_x,mpi_info.mpirank_y,mpi_info.w_rank,mpi_info.e_rank,mpi_info.s_rank,mpi_info.n_rank);
    grid_size = (nx+2) * (ny+2);
    dx = length_x / nx;
    dy = length_y / ny;

    grid_size_mpi = (mpi_info.nx_mpi + 2) * (mpi_info.ny_mpi + 2);

    pos_x = (double*)malloc((mpi_info.nx_mpi+2)*sizeof(double));
    pos_y = (double*)malloc((mpi_info.ny_mpi+2)*sizeof(double));
    u_exact = (double*)malloc(grid_size_mpi*sizeof(double));
    u_solve = (double*)malloc(grid_size_mpi*sizeof(double));
    rhs = (double*)malloc(grid_size_mpi*sizeof(double));

    for(i=0;i<=mpi_info.nx_mpi+1;i++)
    {
        pos_x[i]=(i-0.5 + mpi_info.mpirank_x*mpi_info.nx_mpi)*dx;
    }
    for(j=0;j<=mpi_info.ny_mpi+1;j++)
    {
        pos_y[j]=(j-0.5 + mpi_info.mpirank_y*mpi_info.ny_mpi)*dy;
    }

    for(i=1;i<=mpi_info.nx_mpi;i++)
    {
        x_val = pos_x[i]*(1.0-pos_x[i]);
        for(j=1;j<=mpi_info.ny_mpi;j++)
        {
            y_val = cos(2.0*PI*pos_y[j]);
            u_exact[i*(mpi_info.ny_mpi+2)+j] = x_val * y_val;
            u_solve[i*(mpi_info.ny_mpi+2)+j] = 0.0; 
            rhs[i*(mpi_info.ny_mpi+2)+j] = -2.0*y_val-4.0*PI*PI*x_val*y_val;
        }
    }

    RB_gauss_seidel(dx,dy,mpi_info.nx_mpi,mpi_info.ny_mpi,u_solve,rhs,maxiter,tol,&mpi_info);

    if(mpi_info.myrank == 0)
    {
        FILE *out;

        out = fopen("solution.dat","w+");
        for(i=1;i<=mpi_info.nx_mpi;i++)
        {
            for(j=1;j<=mpi_info.ny_mpi;j++)
            {
                fprintf(out,"%f %f %18.12f %18.12f %18.12f\n",pos_x[i],pos_y[j],rhs[i*(mpi_info.ny_mpi+2)+j],u_exact[i*(mpi_info.ny_mpi+2)+j],u_solve[i*(mpi_info.ny_mpi+2)+j]);
            }
        }
        fclose(out);
    }

    free(pos_x);
    free(pos_y);
    free(u_exact);
    free(u_solve);
    free(rhs);

    MPI_Finalize();

    return 0;
}

