#include "euler.h"
#include "force.h"
#include <stdio.h>
#include <vector>
#include <iostream>

#define cudaH2D cudaMemcpyHostToDevice
#define cudaD2H cudaMemcpyDeviceToHost
#define cudaD2D cudaMemcpyDeviceToDevice

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

void snapshot(vector<triple> save, double *x, double *y, double *s)
{
	double *first = new double[size];
	double *second = new double[size];
	double *third = new double[2];

	checkError(cudaMemcpy(first,x,sizeof(double)*size,cudaD2H));
	checkError(cudaMemcpy(second,y,sizeof(double)*size,cudaD2H));
	checkError(cudaMemcpy(third,s,sizeof(double)*2,cudaD2H));

	triple temp_trip(first,second,third);
	save.push_back(temp_trip);
}

__global__
void vector_sum(double * const out_x,
								double * const out_y,
								double * const out_s,
								const double * const a_x,
								const double * const a_y,
								const double * const a_s,
								const double * const b_x,
								const double * const b_y,
								const double * const b_s)
{
	int index = threadIdx.x;
	out_x[index] = a_x[index] + b_x[index];
	out_y[index] = a_y[index] + b_y[index];
	out_s[index] = a_s[index] + b_s[index];
}

__global__
void update(double * const update_x,
						double * const update_y,
						double * const update_s,
						double * const force_x,
						double * const force_y,
						double * const force_s,
						double h)
{
	int index = threadIdx.x;
	double x = update_x[index];
	double y = update_y[index];
	double s = update_s[index];

	update_x[index] = x + h*force_x[index];
	update_y[index] = y + h*force_y[index];
	update_s[index] = s + h*force_s[index];
}

std::vector<triple>
eulers_method(const double * const x,
							const double * const y,
							const double * const s,
							const double * const delta,
							double t_start,
							double t_end,
							double h,
							int save,
							parameter p)
{
	std::vector<triple> store;

	double *force_x, *force_y, *force_s;
	double *update_x, *update_y, *update_s;
	int size = p.n * p.m;

	checkError(cudaMalloc(&force_x,sizeof(double)*size));
	checkError(cudaMalloc(&force_y,sizeof(double)*size));
	checkError(cudaMalloc(&force_s,sizeof(double)*2));
	checkError(cudaMalloc(&update_x,sizeof(double)*size));
	checkError(cudaMalloc(&update_y,sizeof(double)*size));
	checkError(cudaMalloc(&update_s,sizeof(double)*2));

	checkError(cudaMemcpy(update_x,x,sizeof(double)*size,cudaD2D));
	checkError(cudaMemcpy(update_y,y,sizeof(double)*size,cudaD2D));
	checkError(cudaMemcpy(update_s,s,sizeof(double)*2,cudaD2D));

	double t = t_start;

	int total_points = (int)((t_end - t_start)/h + 1);
	int sampling_rate = total_points / save;

	if ( sampling_rate <= 0 )
	{
		sampling_rate = 1;
	}
	int count = 0;

	cudaStream_t s1, s2;
	checkError(cudaStreamCreate(&s1));
	checkError(cudaStreamCreate(&s2));

	while(t <= t_end)
	{
		if(count % sampling_rate == 0)
		{
			snapshot(store,update_x,update_y,update_s);

			printf("current time: %f, final time: %f",t,t_end);
			std::cout << std::endl;
		}

		force(force_x,force_y,force_s,update_x,update_y,update_s,delta,p);
		update<<<1,size>>>(update_x,update_y,update_s,force_x,force_y,force_s,h);

		++count;
		t+=h;
	}

	snapshot(store,update_x,update_y,update_s);

	checkError(cudaStreamDestroy(s1));
	checkError(cudaStreamDestroy(s2));

	checkError(cudaFree(force_x));
	checkError(cudaFree(force_y));
	checkError(cudaFree(force_s));
	checkError(cudaFree(update_x));
	checkError(cudaFree(update_y));
	checkError(cudaFree(update_s));

	return store;
}


std::vector<triple>
euler_heun_adaptive(const double * const x,
										const double * const y,
										const double * const s,
										const double * const delta,
										double t_start,
										double t_end,
										double h,
										int save,
										double tolerance,
										parameter p)
{
	std::vector<triple> store;

	double *temp_x, *temp_y, *temp_s;
	double *update_x, *update_y, *update_s;
	int size = p.n * p.m;

	checkError(cudaMalloc(&temp_x,sizeof(double)*size));
	checkError(cudaMalloc(&temp_y,sizeof(double)*size));
	checkError(cudaMalloc(&temp_s,sizeof(double)*2));
	checkError(cudaMalloc(&update_x,sizeof(double)*size));
	checkError(cudaMalloc(&update_y,sizeof(double)*size));
	checkError(cudaMalloc(&update_s,sizeof(double)*2));

	checkError(cudaMemcpy(update_x,x,sizeof(double)*size,cudaD2D));
	checkError(cudaMemcpy(update_y,y,sizeof(double)*size,cudaD2D));
	checkError(cudaMemcpy(update_s,s,sizeof(double)*2,cudaD2D));

	double *k1_x, *k1_y, *k1_s, *k2_x, *k2_y, *k2_s;
	checkError(cudaMalloc(&k1_x,sizeof(double)*size));
	checkError(cudaMalloc(&k1_y,sizeof(double)*size));
	checkError(cudaMalloc(&k1_s,sizeof(double)*size));
	checkError(cudaMalloc(&k2_x,sizeof(double)*size));
	checkError(cudaMalloc(&k2_y,sizeof(double)*size));
	checkError(cudaMalloc(&k2_s,sizeof(double)*size));

	double t = t_start;

	int total_points = (int)((t_end - t_start)/h + 1);
	int sampling_rate = total_points / save;

	if ( sampling_rate <= 0 )
	{
		sampling_rate = 1;
	}
	int count = 0;

	double error = 100;
	while(t <= t_end)
	{
		if(count % sampling_rate == 0)
		{
			snapshot(store,update_x,update_y,update_s);

			printf("current time: %f, final time: %f",t,t_end);
			std::cout << std::endl;
		}

		while(error > tolerance)
		{
			force(k1_x,k1_y,k1_s,update_x,update_y,update_s,delta,p);
			update<<<1,size>>>(temp_x,temp_y,temp_s,k1_x,k1_y,k1_s,h);
			force(k2_s,k2_y,k2_s,temp_x,temp_y,temp_s,delta,p);
			vector_sum<<<1,size>>>(temp_x,temp_y,temp_s,k1_x,k1_y,k1_s,k2_x,k2_y,k2_s);
			update<<<1,size>>>(update_x,update_y,update_s,temp_x,temp_y,temp_s,h/2);
		}

		++count;
		t+=h;
	}

	snapshot(store,update_x,update_y,update_s);

	checkError(cudaFree(temp_x));
	checkError(cudaFree(temp_y));
	checkError(cudaFree(temp_s));
	checkError(cudaFree(update_x));
	checkError(cudaFree(update_y));
	checkError(cudaFree(update_s));
	checkError(cudaFree(k1_x));
	checkError(cudaFree(k1_y));
	checkError(cudaFree(k1_s));
	checkError(cudaFree(k2_x));
	checkError(cudaFree(k2_y));
	checkError(cudaFree(k2_s));

	return store;
}
