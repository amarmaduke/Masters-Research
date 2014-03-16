#include "jacobi.h"

__device__ void position(int& j, int& i, const int& x, const parameter& p)
{
	j = x / p.n;
	i = (x % p.n);
}

__global__ void naive_jacobi(double* matrix,
							const size_t pitch,
							const double* x,
							const double* y,
							const double* s,
							const parameter p)
{
	int r = threadIdx.x, c = threadIdx.y;
	int j, i, h, k;
	int outptr = r*(pitch/sizeof(value_type)) + c;
	//printf("r: %d, c: %d, pitch: %d, out: %d, outptr: %d\n",r,c,pitch,(r*pitch+c),outptr);
	position(j,i,r,p);
	position(h,k,c,p);
	printf("r: %d, c: %d, outptr: %d, j: %d, i: %d\n",r,c,outptr,j,i);
	matrix[outptr] = 1;
	if(r != c)
		matrix[outptr] = 2;
}

void jacobi_functor::operator() (const vector_type& x, matrix_type& J, const value_type& t, vector_type& dfdt)
{
	int rows = J.size1(), cols = J.size2();
	int n = state.n, m = state.m;
	int fiber_size = n*m;
	dim3 grid_size(1,1,1), block_size(rows,cols,1);
	double* h_matrix = &(J.data()[0]);
	const double* h_x = &(x.data()[0]);
	size_t h_pitch = cols*sizeof(value_type);
	double* d_matrix, * d_x;
	size_t d_pitch;
	cudaMallocPitch(&d_matrix,&d_pitch,cols*sizeof(value_type),rows);
	cudaMalloc(&d_x,x.size()*sizeof(double));
	cudaMemcpy(d_x,h_x,x.size()*sizeof(double),cudaMemcpyHostToDevice);

	naive_jacobi<<<grid_size,block_size>>>
		(d_matrix,d_pitch,d_x,d_x+fiber_size,d_x+2*fiber_size,state);


	cudaMemcpy2D(	h_matrix, h_pitch,
								d_matrix, d_pitch,
								cols*sizeof(value_type), rows, cudaMemcpyDeviceToHost);

}


