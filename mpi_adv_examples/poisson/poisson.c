#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void RB_gauss_seidel(double dx, double dy, int grid_x, int grid_y, double *u_solve, double *rhs, int maxiteration, double tolerance)
{
	double *u_new;
	double error_sum, u_sum;
	double ae,aw,as,an,ap;
	double dxsqi,dysqi;
	double uij_old,uij_new;
	int i,j,iter,js,walk;

	dxsqi=1.0/dx/dx;
	dysqi=1.0/dy/dy;

	for(iter=0;iter<maxiteration;iter++)
	{
		error_sum = 0.0;
		u_sum = 0.0;

		for(walk=1;walk<3;walk++)
		{
			for(i=1;i<=grid_x;i++)
			{
				u_solve[i*(grid_y+2)]=u_solve[i*(grid_y+2)+grid_y];
				u_solve[i*(grid_y+2)+grid_y+1]=u_solve[i*(grid_y+2)+1];
			}
			for(j=1;j<=grid_y;j++)
			{
				u_solve[0*(grid_y+2)+j]=-u_solve[1*(grid_y+2)+j];
				u_solve[(grid_x+1)*(grid_y+2)+j]=-u_solve[grid_x*(grid_y+2)+j];
			}
			js = 3 - walk;
			for(i=1;i<=grid_x;i++)
			{
				js=3-js;
				for(j=js;j<=grid_y;j+=2)
				{
					aw=dxsqi*u_solve[(i-1)*(grid_y+2)+j];
					ae=dxsqi*u_solve[(i+1)*(grid_y+2)+j];
					as=dysqi*u_solve[i*(grid_y+2)+j-1];
					an=dysqi*u_solve[i*(grid_y+2)+j+1];
	
					ap = 2.0*(dxsqi+dysqi);
	
					uij_old = u_solve[i*(grid_y+2)+j];
					uij_new = (aw+ae+as+an-rhs[i*(grid_y+2)+j]) /ap;
	
					error_sum += fabs(uij_new-uij_old);
					u_sum += fabs(uij_new);
	
					u_solve[i*(grid_y+2)+j] = uij_new;
				}
			}
		}
		if(iter%1==0) printf("%5d th iteration: error=%e\n",iter,error_sum/u_sum);
		if(error_sum/u_sum<tolerance) break;
	}
	printf("Iteration ends in %d th step\n",iter);
}

int main(int argc, char **argv)
{
	int grid_x = 256, grid_y = 256;
        int maxiter = 100000;
        double tol = 1.0e-8;
	double length_x = 1.0, length_y = 1.0;
	double PI = atan(1.0)*4.0;

	int grid_size;
	double dx, dy, x_val, y_val;
	double *pos_x, *pos_y; double *u_exact, *u_solve, *rhs; 
	int i,j;

	grid_size = (grid_x+2) * (grid_y+2);
	dx = length_x / grid_x;
	dy = length_y / grid_y;

	pos_x = (double*)malloc((grid_x+2)*sizeof(double));
	pos_y = (double*)malloc((grid_y+2)*sizeof(double));
	u_exact = (double*)malloc(grid_size*sizeof(double));
	u_solve = (double*)malloc(grid_size*sizeof(double));
	rhs = (double*)malloc(grid_size*sizeof(double));

	for(i=0;i<=grid_x+1;i++)
	{
		pos_x[i]=(i-0.5)*dx;
	}
	for(j=0;j<=grid_y+1;j++)
	{
		pos_y[j]=(j-0.5)*dy;
	}

	for(i=1;i<=grid_x;i++)
	{
		x_val = pos_x[i]*(1.0-pos_x[i]);
		for(j=1;j<=grid_y;j++)
		{
			y_val = cos(2.0*PI*pos_y[j]);
			u_exact[i*(grid_y+2)+j] = x_val * y_val;
			u_solve[i*(grid_y+2)+j] = 0.0; 
			rhs[i*(grid_y+2)+j] = -2.0*y_val-4.0*PI*PI*x_val*y_val;
		}
	}

	RB_gauss_seidel(dx,dy,grid_x,grid_y,u_solve,rhs,maxiter,tol);

	FILE *out;

	out = fopen("solution.dat","w+");
	for(i=1;i<=grid_x;i++)
	{
		for(j=1;j<=grid_y;j++)
		{
			fprintf(out,"%f %f %18.12f %18.12f %18.12f\n",pos_x[i],pos_y[j],rhs[i*(grid_y+2)+j],u_exact[i*(grid_y+2)+j],u_solve[i*(grid_y+2)+j]);
		}
	}
	fclose(out);

	free(pos_x);
	free(pos_y);
	free(u_exact);
	free(u_solve);
	free(rhs);


	return 0;
}

