#include "force.h"
#include "parameter.h"
#include <iostream>
#include <stdio.h>

#define _DEBUG_

#ifdef _DEBUG_
	#define debug(x) x
#else
	#define debug(x)
#endif

#ifdef _ERROR_
#define checkError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n",
		cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}
#else
#define checkError(ans) ans
#endif


__global__
void compute_lengths(	double * const lens, 
											const double * const in_x,
											const double * const in_y,
											const parameter * const p)
{
	int j = threadIdx.x;
	int i = threadIdx.y;
	int n = p->n;
	const double * const delta = p->delta;

	int index = j*n + i;

	double x = in_x[index];
	double y = in_y[index];

	double xp = i == 0 ? delta[j] : in_x[index-1];
	double yp = i == 0 ? 0 : in_y[index-1];

	double l = (xp-x)*(xp-x)+(yp-y)*(yp-y);

	lens[index] = l;

	debug(printf(
		"block(%d,%d,%d), thread(%d,%d,%d), index:%d\n"
		"x: %f, y: %f\n"
		"xp: %f, yp: %f\n"
		"l: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,
		threadIdx.x,threadIdx.y,threadIdx.z,index,x,y,xp,yp,l
	);)
}

__global__
void compute_fiber_dependent(	double * const out_x,
															double * const out_y,
															const double * const in_x,
															const double * const in_y,
															const double * const lens,
															const parameter * const p)
{
	int j = threadIdx.x; // What fiber we're on
	int i = threadIdx.y; // What particle it is
	int n = p->n;
	double beta = p->beta;
	double len = p->len;
	double gamma = p->gamma;
	double epsilon = p->epsilon;
	double sigma = p->sigma;
	const double* const delta = p->delta;

	int index = j*n + i;

	double xp = i == 0 ? delta[j] : in_x[index-1];
	double xpp = i == 0 || i == 1 ? delta[j] : in_x[index-2];

	double yp = i == 0 ? 0 : in_y[index-1];
	double ypp = 	i == 0 || i == 1 ?
							( i == 0 ? -len : 0) : in_y[index-2];

	double xn = (index % n) + 1 < n ? in_x[index+1] : NAN;
	double xnn = (index % n) + 2 < n ? in_x[index+2] : NAN;
	double yn = (index % n) + 1 < n ? in_y[index+1] : NAN;
	double ynn = (index % n) + 2 < n ? in_y[index+2] : NAN;

	double l = lens[index];
	double x = in_x[index];
	double y = in_y[index];

	double lp = i == 0 ? sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
							: lens[index-1];
	double ln = (index % n) + 1 < n ? lens[index+1] : NAN;
	double lnn = (index % n) + 2 < n ? lens[index+2] : NAN;

	// Compute Bending Force

	double xd_f = (xn - x)*(xnn - xn);
	double yd_f = (yn - y)*(ynn - yn);
	double xd_c = (x - xp)*(xn - x);
	double yd_c = (y - yp)*(yn - y);
	double xd_b = (xp - xpp)*(x - xp);
	double yd_b = (yp - ypp)*(y - yp);

	double product_f = xd_f + yd_f;
	double product_c = xd_c + yd_c;
	double product_b = xd_b + yd_b;

	double forward_1 = lnn / (ln * product_f);
	double forward_2 = ln*lnn/(product_f*product_f);
	double center_1 = ln/(l*product_c);
	double center_2 = l/(ln*product_c);
	double center_3 = l*ln/(product_c*product_c);
	double backward_1 = lp/(l*product_b);
	double backward_2 = lp*l/(product_b*product_b);

	double forward_1x_t = forward_1*(x-xn);
	double forward_2x_t = forward_2*(xnn-xn);
	double center_1x_t = center_1*(x-xp);
	double center_2x_t = center_2*(x-xn);
	double center_3x_t = center_3*(xn-2*x+xp);
	double backward_1x_t = backward_1*(x-xp);
	double backward_2x_t = backward_2*(xpp-xp);

	double forward_1y_t = forward_1*(y-yn);
	double forward_2y_t = forward_2*(ynn-yn);
	double center_1y_t = center_1*(y-yp);
	double center_2y_t = center_2*(y-yn);
	double center_3y_t = center_3*(yn-2*y+yp);
	double backward_1y_t = backward_1*(y-yp);
	double backward_2y_t = backward_2*(ypp-yp);

	double forward_1x = isnan(forward_1x_t) ? 0 : forward_1x_t;
	double forward_2x = isnan(forward_2x_t) ? 0 : forward_2x_t;
	double center_1x = isnan(center_1x_t) ? 0 : center_1x_t;
	double center_2x = isnan(center_2x_t) ? 0 : center_2x_t;
	double center_3x = isnan(center_3x_t) ? 0 : center_3x_t;
	double backward_1x = isnan(backward_1x_t) ? 0 : backward_1x_t;
	double backward_2x = isnan(backward_2x_t) ? 0 : backward_2x_t;

	double forward_1y = isnan(forward_1y_t) ? 0 : forward_1y_t;
	double forward_2y = isnan(forward_2y_t) ? 0 : forward_2y_t;
	double center_1y = isnan(center_1y_t) ? 0 : center_1y_t;
	double center_2y = isnan(center_2y_t) ? 0 : center_2y_t;
	double center_3y = isnan(center_3y_t) ? 0 : center_3y_t;
	double backward_1y = isnan(backward_1y_t) ? 0 : backward_1y_t;
	double backward_2y = isnan(backward_2y_t) ? 0 : backward_2y_t;

	double forward_x = forward_1x + forward_2x;
	double center_x = center_1x + center_2x + center_3x;
	double backward_x = backward_1x + backward_2x;

	double forward_y = forward_1y + forward_2y;
	double center_y = center_1y + center_2y + center_3y;
	double backward_y = backward_1y + backward_2y;

	double bending_x = beta*(forward_x + center_x + backward_x);
	double bending_y = beta*(forward_y + center_y + backward_y);

	// Compute Extensible Spring Force

	double e_forward = (ln - len)/ln;
	double e_backward = (l - len)/l;

	double e_forward_x_t = e_forward*2*(x-xn);
	double e_backward_x_t = e_backward*2*(x-xp);
	double e_forward_y_t = e_forward*2*(y-yn);
	double e_backward_y_t = e_backward*2*(y-yp);

	double e_forward_x = isnan(e_forward_x_t) ? 0 : e_forward_x_t;
	double e_backward_x = isnan(e_backward_x_t) ? 0 : e_backward_x_t;
	double e_forward_y = isnan(e_forward_y_t) ? 0 : e_forward_y_t;
	double e_backward_y = isnan(e_backward_y_t) ? 0 : e_backward_y_t;

	double extensible_x = gamma*(e_forward_x + e_backward_x);
	double extensible_y = gamma*(e_forward_y + e_backward_y);

	// Compute particle-substrate vdW

	double p1 = sigma / y;
	double p2 = p1*p1;
	double p4 = p2*p2;
	double p5 = p4*p1;
	double p11 = p5*p5*p1;

	double vdW_y = -(3.14*epsilon)*(2*p11-4*p5);

	// Compute Total Force

	double total_force_x = -(bending_x + extensible_x);
	double total_force_y = -(bending_y + extensible_y + vdW_y);

	///*
	debug(printf(
		"compute fiber dependent force:\n"
		"Block index: (%d,%d,%d), Thread index: (%d,%d,%d)\n"
		"xpp: %.5f, xp: %.5f, x: %.5f, xn: %.5f, xnn: %.5f\n"
		"ypp: %.5f, yp: %.5f, y: %.5f, yn: %.5f, ynn: %.5f\n"
		"lp: %.5f, l: %.5f, ln: %.5f, lnn: %.5f\n"
		"bending_x: %f, bending_y: %f\n"
		"extensible_x: %f, extensible_y: %f\n"
		"vdW_y: %f\n"
		"total_x: %f, total_y: %f\n",
		blockIdx.x,blockIdx.y,blockIdx.z,
		threadIdx.x,threadIdx.y,threadIdx.z,
		xpp,xp,x,xn,xnn,ypp,yp,y,yn,ynn,lp,l,ln,lnn,
		bending_x,bending_y,extensible_x,extensible_y,
		vdW_y,total_force_x,total_force_y
	);) //*/

	out_x[index] = total_force_x;
	out_y[index] = total_force_y;
}



__global__
void compute_n_body_vdw(double * const out_x,
												double * const out_y,
												double * const out_s,
												const double * const in_x,
												const double * const in_y,
												const double * const in_s,
												const parameter * const p)
{
	int j = threadIdx.x; // What fiber we're on
	int i = threadIdx.y; // What particle it is

	double n = p->n;
	double m = p->m;
	double epsilon = p->epsilon;
	double sigma = p->sigma;
	double sub_h = p->sub_h;
	double sub_count = p->sub_count;

	int index = j*n + i;

	double x = in_x[index];
	double y = in_y[index];

	double vdW_x = 0, vdW_y = 0;
	for(int j_ = 0; j_ < m; ++j_)
	{
		for(int i_ = 0; i_ < n; ++i_)
		{
			if(j == j_ && (i_ == i || i_ + 1 == i || i_ - 1 == i))
				continue;
			int k = j_*n + i_;

			double x_ = in_x[k];
			double y_ = in_y[k];

			double xps = x - x_;
			double yps = y - y_;
			double dist = sqrt(xps*xps + yps*yps);

			double temp_x = xps/dist;
			double temp_y = yps/dist;

			double p1 = sigma / dist;
			double p2 = p1*p1;
			double p4 = p2*p2;
			double p7 = p4*p2*p1;
			double p8 = p7*p1;
			double p13 = p8*p4*p1;
			double LJval = -(12*epsilon/sigma)*(p13-p7);

			vdW_x = vdW_x - LJval*temp_x;
			vdW_y = vdW_y - LJval*temp_y;

			///*
			debug(printf(
				"Index: %d, Index:o: %d\n"
				"vdW_x: %f, vdW_y: %f\n",
				index,k,LJval*temp_x,LJval*temp_y
			);) //*/
		}
	}

	double s_x = in_s[0];
	double s_y = in_s[1];
	double s_vdW_x = 0, s_vdW_y = 0;

	for(int k = 0; k < sub_count; ++k)
	{
		double x_ = s_x +k*sub_h;
		double y_ = s_y;

		double xps = x - x_;
		double yps = y - y_;
		double dist = sqrt(xps*xps + yps*yps);

		double temp_x = xps/dist;
		double temp_y = yps/dist;

		double p1 = sigma / dist;
		double p2 = p1*p1;
		double p4 = p2*p2;
		double p7 = p4*p2*p1;
		double p8 = p7*p1;
		double p13 = p8*p4*p1;
		double LJval = -(12*epsilon/sigma)*(p13-p7);

		vdW_x = vdW_x - LJval*temp_x;
		vdW_y = vdW_y - LJval*temp_y;
		s_vdW_x = s_vdW_x + LJval*temp_x;
		s_vdW_y = s_vdW_y + LJval*temp_y;
	}

	out_s[0] = out_s[0] + s_vdW_x;
	out_s[1] = out_s[1] + s_vdW_y;

	out_x[index] = vdW_x;
	out_y[index] = vdW_y;
}

__global__
void combine(	double * const out_x,
							double * const out_y,
							const double * const f_1x,
							const double * const f_1y,
							const double * const f_2x,
							const double * const f_2y,
							const parameter * const p)
{
	int j = threadIdx.x; // What fiber we're on
	int i = threadIdx.y; // What particle it is

	double n = p->n;
	int index = j*n + i;

	out_x[index] = f_1x[index] + f_2x[index];
	out_y[index] = f_1y[index] + f_2y[index];

	debug(printf(
		"Index: %d\n"
		"force_x: %f\n"
		"force_y: %f\n",
		index,out_x[index],out_y[index]
	);)
}


__global__
void print_all(double * out_x,
							 double * out_y,
							 double * out_s,
							 const double * in_x,
							 const double * in_y,
							 const double * in_s,
							 double * lens,
							 double * f_1x,
							 double * f_1y,
							 double * f_2x,
							 double * f_2y,
							 const parameter * p)
{
	int size = p->m * p->n;

	printf("out_x: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",out_x[i]);
	} printf("\n");

	printf("out_y: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",out_y[i]);
	} printf("\n");

	printf("out_s: ");
	for(int i = 0; i < 2; ++i)
	{
		printf("%f, ",out_x[i]);
	} printf("\n");

	printf("in_x: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",in_x[i]);
	} printf("\n");

	printf("in_y: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",in_y[i]);
	} printf("\n");

	printf("in_s: ");
	for(int i = 0; i < 2; ++i)
	{
		printf("%f, ",in_s[i]);
	} printf("\n");

	printf("lens: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",lens[i]);
	} printf("\n");

	printf("f_1x: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",f_1x[i]);
	} printf("\n");

	printf("f_1y: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",f_1y[i]);
	} printf("\n");

	printf("f_2x: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",f_2x[i]);
	} printf("\n");

	printf("f_2y: ");
	for(int i = 0; i < size; ++i)
	{
		printf("%f, ",f_2y[i]);
	} printf("\n");

}


void force(	double * const out_x,
						double * const out_y,
						double * const out_s,
						const double * const in_x,
						const double * const in_y,
						const double * const in_s,
						const parameter h_p,
						const parameter * const d_p)
{
	dim3 grid(1,1,1), blocks(h_p.m,h_p.n,1);
	int size = h_p.n * h_p.m;

	printf("Hello World\n");

	double *lens, *f_1x, *f_1y, *f_2x, *f_2y;
	checkError(cudaMalloc(&lens,sizeof(double)*size));
	checkError(cudaMalloc(&f_1x,sizeof(double)*size));
	checkError(cudaMalloc(&f_1y,sizeof(double)*size));
	checkError(cudaMalloc(&f_2x,sizeof(double)*size));
	checkError(cudaMalloc(&f_2y,sizeof(double)*size));

	checkError(cudaMemset(out_s,0,sizeof(double)*2));
	checkError(cudaMemset(lens,0,sizeof(double)*size));
	checkError(cudaMemset(out_x,0,sizeof(double)*size));
	checkError(cudaMemset(out_y,0,sizeof(double)*size));
	checkError(cudaMemset(f_1x,0,sizeof(double)*size));
	checkError(cudaMemset(f_1y,0,sizeof(double)*size));
	checkError(cudaMemset(f_2x,0,sizeof(double)*size));
	checkError(cudaMemset(f_2y,0,sizeof(double)*size));


	print_all<<<1,1>>>(out_x,out_y,out_s,in_x,in_y,in_s,lens,f_1x,f_1y,f_2x,f_2y,d_p);

	//printf("Hello World\n");
	int i = 0;


	compute_lengths<<<grid,blocks>>>(lens,in_x,in_y,d_p);
	compute_fiber_dependent<<<grid,blocks>>>(f_1x,f_1y,in_x,in_y,lens,d_p);
	compute_n_body_vdw<<<grid,blocks>>>(f_2x,f_2y,out_s,in_x,in_y,in_s,d_p);
	combine<<<grid,blocks>>>(out_x,out_y,f_1x,f_1y,f_2x,f_2y,d_p);

	print_all<<<1,1>>>(out_x,out_y,out_s,in_x,in_y,in_s,lens,f_1x,f_1y,f_2x,f_2y,d_p);

	cudaDeviceSynchronize();
	scanf("%d",&i);

	checkError(cudaFree(lens));
	checkError(cudaFree(f_1x));
	checkError(cudaFree(f_1y));
	checkError(cudaFree(f_2x));
	checkError(cudaFree(f_2y));
}

