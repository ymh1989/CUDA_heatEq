#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "utils.h"
#include "dev_matrix.h"

// DIVIDE_INTO(x/y) for integers, used to determine # of blocks/warps etc.
#define DIVIDE_INTO(x,y) (((x) + (y) - 1)/(y))
// I2D to index into a linear memory space from a 2D array index pair
#define I2D(Nx, i, j) ((i) + (Nx)*(j))

// Block size in the i and j directions
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16

// kernel to update temperatures - CPU version
void heat2d_cpu(int Nx, int Ny, double alp, double *in, double *out) {
	int i, j, P, W, E, S, N;
	double d2tdx2, d2tdy2;
	// loop over all points in domain (not boundary points)
	for (j = 1; j < Ny-1; j++) {
		for (i = 1; i < Nx-1; i++) {
			// find indices into linear memory for central point and neighbours
			P = I2D(Nx, i, j);
			W = I2D(Nx, i-1, j); E = I2D(Nx, i+1, j);
			S = I2D(Nx, i, j-1); N = I2D(Nx, i, j+1);

			d2tdx2 = in[W] - 2.0*in[P] + in[E];
			d2tdy2 = in[N] - 2.0*in[P] + in[S];

			out[P] = in[P] + alp*(d2tdx2 + d2tdy2);
		}
	}
}

// kernel to update temperatures - GPU version (not using shared mem)
__global__ void heat2d_gpu(int Nx, int Ny, double alp, double *in, double *out) 
{
	int i, j, P, W, E, S, N;
	double d2tdx2, d2tdy2;
	// find i and j indices of this thread
	i = blockIdx.x*(BLOCK_SIZE_X) + threadIdx.x;
	j = blockIdx.y*(BLOCK_SIZE_Y) + threadIdx.y;

	// find indices into linear memory 
	P = I2D(Nx, i, j);
	W = I2D(Nx, i-1, j); E = I2D(Nx, i+1, j);
	S = I2D(Nx, i, j-1); N = I2D(Nx, i, j+1);

	// check that thread is within domain (not on boundary or outside domain)
	if (i > 0 && i < Nx-1 && j > 0 && j < Ny-1) {
		d2tdx2 = in[W] - 2.0*in[P] + in[E];
		d2tdy2 = in[N] - 2.0*in[P] + in[S];

		out[P] = in[P] + alp*(d2tdx2 + d2tdy2);
	}
}

void heat2d_exc(double * out, double * in1, double * in2,
	const double T, const int Nx, const int Ny)
{		
	int i, j;

	for (j = 1; j < Ny-1; j++) {
		for (i = 1; i < Nx-1; i++) {
			out[I2D(Nx, i, j)] = sin(M_PI*in1[i])*sin(M_PI*in2[j]) * exp(-2.0*M_PI*M_PI*T);
		}
	}
}


int main() 
{
	int Nx, Ny, Nt;
	double alpha, *x, *y, *u_h, *oldu_h, *tmp_h, *exc;
	double * u_d;
	int i, j, iter;
	double h, T;
	double errGPU = 0.0, errCPU = 0.0;
	dim3 numBlocks, threadsPerBlock;
	double clock_h, clock_d;
	FILE *fp;

	// domain size and number of timesteps (iterations)
	Nx = 1024;
	Ny = Nx;
	Nt = 100;
	alpha = 0.25;
	h = 1.0 / (double)(Nx-1);
	T = Nt*alpha*h*h;

	// allocate temperature array on host
	// x = (double *)malloc(sizeof(double)*Nx);
	x = dvector(Nx); y = dvector(Ny);
	u_h = dvector(Nx*Ny); oldu_h = dvector(Nx*Ny);
	exc = dvector(Nx*Ny);
	u_d = dvector(Nx*Ny);

	for (j = 0; j < Nx-1; j++) x[j] = 0.0 + (j*h);
	for (i = 0; i < Ny-1; i++) y[i] = 0.0 + (i*h);

	zero_matrix(u_h, Nx, Ny);
	zero_matrix(oldu_h, Nx, Ny);
	zero_matrix(exc, Nx, Ny);

	// initial
	initialize(u_h, x, y, Nx, Ny);
	initialize(oldu_h, x, y, Nx, Ny);

	// allocate temperature arrays on device
	dev_matrix<double> ud(Nx, Ny); ud.set(u_h, Nx, Ny);
	dev_matrix<double> oldud(Nx, Ny); oldud.set(u_h, Nx, Ny);
	dev_matrix<double> tmp_d(Nx, Ny);

	// set threads and blocks
	numBlocks = dim3(iDivUp(Nx,BLOCK_SIZE_X), iDivUp(Ny,BLOCK_SIZE_Y));
	threadsPerBlock = dim3(BLOCK_SIZE_X, BLOCK_SIZE_Y);

	// cpu loop
	printf("CPU start!\n");
	clock_h = double(clock()) / CLOCKS_PER_SEC;
	for (iter = 0; iter < Nt; iter++) {
		heat2d_cpu(Nx, Ny, alpha, oldu_h, u_h);
        tmp_h = u_h;
        u_h = oldu_h;
        oldu_h = tmp_h;
	}
	clock_h = double(clock()) / CLOCKS_PER_SEC - clock_h;
	printf("CPU end!\n");

	// gpu loop
	printf("GPU start!\n");
	clock_d = double(clock()) / CLOCKS_PER_SEC;
	for (iter = 0; iter < Nt; iter++) {
		heat2d_gpu<<<numBlocks, threadsPerBlock>>>(Nx, Ny, alpha, oldud.getData(), ud.getData());
        tmp_d = ud;
        ud = oldud;
        oldud = tmp_d;
	} 
	cudaThreadSynchronize();
	clock_d = double(clock()) / CLOCKS_PER_SEC - clock_d;
	printf("GPU end!\n");

	// copy temperature array from device to host
	oldud.get(&u_d[0], Nx, Ny);

	// Exact value
	heat2d_exc(exc, x, y, T, Nx, Ny);

	//RMSE
	for (i = 0; i < Nx*Ny; i++) {			
		errGPU = errGPU + ( (u_d[i]-exc[i])*(u_d[i]-exc[i]) );
		errCPU = errCPU + ( (oldu_h[i]-exc[i])*(oldu_h[i]-exc[i]) );
	}
	errGPU = sqrt(errGPU / (Nx*Ny));
	errCPU = sqrt(errCPU / (Nx*Ny));
	printf("\n");
	printf("RMSE (CPU) : %.12f\n", errCPU);
	printf("RMSE (GPU) : %.12f\n", errGPU);
	printf("CPU time = %.3fms\n",clock_h*1e3);
	printf("GPU time = %.3fms\n",clock_d*1e3);
	printf("CPU time / GPU time : %.2f\n", clock_h/clock_d);

	printf("\n");
	printf("Printing...\n");
	fp = fopen("host_out.dat", "w");
	print_mat(fp, oldu_h, Nx, Ny);
	
	fp = fopen("dev_out.dat", "w");
	print_mat(fp, u_d, Nx, Ny);
	
	fp = fopen("exact.dat", "w");
	print_mat(fp, exc, Nx, Ny);


	oldud.~dev_matrix(); ud.~dev_matrix();
	free(x); free(y); free(u_h); free(oldu_h);
	free(exc); free(u_d);

	return 0;
}


