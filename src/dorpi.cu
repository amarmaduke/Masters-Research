#include "dorpi.h"
#include "force.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cuda_runtime.h>
#include "cublas_v2.h"

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

void SNAPSHOT(std::vector<triple>& save, double *x, double *y, double *s, int size)
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

template<typename T>
void print_v(int n, T * x)
{
	std::cout << "(";
	for(int i = 0; i < n; ++i)
	{
		std::cout << x[i];
		if(i + 1 != n)
			std::cout << " ";
	}
	std::cout << ")";
}

// count = # of vectors passed to linc
// size = # of elements in each vector
// ... is assumed to have count (constant,vector) pairs
double * linc(cublasHandle_t handle, int count, int size, ...)
{
  va_list ap;
	double ** vectors, * constants;
  int * iter;
	
	vectors = new double*[count];
	constants = new double[count];
	iter = new int[count];

  // Parse variadic input and store copies
  va_start(ap,size);
	for(int i = 0; i < count; ++i)
  {
    // Setup iterator
    iter[i] = i;

    // Grab constant
		double c_temp = va_arg(ap, double);
		constants[i] = c_temp;

    // Grab vector
    checkError(cudaMalloc(&vectors[i], sizeof(double)*size));
		double * v_temp = va_arg(ap, double *);
		
		cublasDcopy(handle, size, v_temp, 1, vectors[i], 1);
	}
  va_end(ap);

  // Compute linear compination
  int N = count, m, j = 0;
  while(N > 1)
  {
		std::cout << "iter: " << std::endl;
		print_v(count,iter);
		std::cout << std::endl;
    m = N % 2 == 0 ? N : N - 1;
    for(int i = 0; i < m; i+=2, ++j)
    {
      double a = constants[iter[i+1]]/constants[iter[i]];
      
			double *x, *y;
			x = new double[size]; y = new double[size];
			cudaMemcpy(x,vectors[iter[i+1]],sizeof(double)*size,cudaD2H);
			cudaMemcpy(y,vectors[iter[i]],sizeof(double)*size,cudaD2H);
			
			std::cout << " Inner linc Loop: " << std::endl;
			std::cout << a;
			print_v(size,x);
			std::cout << " + ";
			print_v(size,y);
			std::cout << std::endl;

			cublasDaxpy(handle, count, &a, vectors[iter[i+1]], 1, vectors[iter[i]], 1);
      iter[j] = iter[i];
    }
    if(N % 2 != 0)
      iter[j] = iter[m];
    N = (N+1)/2;
    j = 0;
		std::cout << " Outer liner loop: " << std::endl;
		std::cout << "N: " << N << std::endl;
  }
	
	double a = constants[iter[0]];
	cublasDscal(handle, count, &a, vectors[0], 1);

  // Grab output and clean
  double * out = vectors[0];
  for(int i = 1; i < count; ++i)
  {
    cudaFree(vectors[i]);
  }
  delete constants;
  delete iter;
	delete vectors;

  return out;
}

void test()
{
	int N = 3;
	double *h_x, *h_y, *h_z, *h_a, *h_b, *h_r;
	double *d_x, *d_y, *d_z, *d_a, *d_b, *d_r;

	h_x = new double[N];
	h_y = new double[N];
	h_z = new double[N];
	h_a = new double[N];
	h_b = new double[N];
	h_r = new double[N];

	for(int i = 0; i < N; ++i)
	{
		h_x[i] = 1;
		h_y[i] = 1;
		h_z[i] = 1;
		h_a[i] = 1;
		h_b[i] = 1;
	}

	cudaMalloc(&d_x, sizeof(double)*N);
	cudaMalloc(&d_y, sizeof(double)*N);
	cudaMalloc(&d_z, sizeof(double)*N);
	cudaMalloc(&d_a, sizeof(double)*N);
	cudaMalloc(&d_b, sizeof(double)*N);
	
	cudaMemcpy(d_x,h_x,sizeof(double)*N,cudaH2D);
	cudaMemcpy(d_y,h_y,sizeof(double)*N,cudaH2D);
	cudaMemcpy(d_z,h_z,sizeof(double)*N,cudaH2D);
	cudaMemcpy(d_a,h_a,sizeof(double)*N,cudaH2D);
	cudaMemcpy(d_b,h_b,sizeof(double)*N,cudaH2D);

	cublasHandle_t handle;
	cublasCreate(&handle);

	std::cout << "1";
	print_v(N,h_x);
	std::cout << std::endl;
	std::cout << "2";
	print_v(N,h_y);
	std::cout << std::endl;
	std::cout << "3";
	print_v(N,h_z);
	std::cout << std::endl;
	std::cout << "2";
	print_v(N,h_a);
	std::cout << std::endl;
	std::cout << "3";
	print_v(N,h_b);
	std::cout << std::endl;


	d_r = linc(handle,5,N,1.0,d_x,2.0,d_y,3.0,d_z,2.0,d_a,3.0,d_b);

	cudaMemcpy(h_r,d_r,sizeof(double)*N,cudaD2H);
	
	for(int i = 0; i < N; ++i)
	{
		std::cout << h_r[i] << " ";
	}
	std::cout << std::endl;
}

__global__
void construct_next_step( double * const out_x,
                          double * const out_y,
                          double * const out_s,
                          const double * const x,
                          const double * const y,
                          const double * const s,
                          const double * const k1_x,
                          const double * const k1_y,
                          const double * const k1_s,
                          const double * const k2_x,
                          const double * const k2_y,
                          const double * const k2_s,
                          const double * const k3_x,
                          const double * const k3_y,
                          const double * const k3_s,
                          const double * const k4_x,
                          const double * const k4_y,
                          const double * const k4_s,
                          const double * const k5_x,
                          const double * const k5_y,
                          const double * const k5_s,
                          const double * const k6_x,
                          const double * const k6_y,
                          const double * const k6_s,
                          const double * const a,
                          const double h)
{
  int index = threadIdx.x;
  out_x[index] = x[index] + h*(a[0]*k1_x[index] + a[1]*k2_x[index]
                 + a[2]*k3_x[index] + a[3]*k4_x[index]
                 + a[4]*k5_x[index] + a[5]*k6_x[index]);
  out_y[index] = y[index] + h*(a[0]*k1_y[index] + a[1]*k2_y[index]
                 + a[2]*k3_y[index] + a[3]*k4_y[index]
                 + a[4]*k5_y[index] + a[5]*k6_y[index]);
  out_s[0] = s[0] + h*(a[0]*k1_s[0] + a[1]*k2_s[0]
                 + a[2]*k3_s[0] + a[3]*k4_s[0]
                 + a[4]*k5_s[0] + a[5]*k6_s[0]);
  out_s[1] = s[1] + h*(a[0]*k1_s[1] + a[1]*k2_s[1]
                 + a[2]*k3_s[1] + a[3]*k4_s[1]
                 + a[4]*k5_s[1] + a[5]*k6_s[1]);
}

__global__
void construct_error_step(double * const out_x,
                          double * const out_y,
                          double * const out_s,
                          const double * const k1_x,
                          const double * const k1_y,
                          const double * const k1_s,
                          const double * const k2_x,
                          const double * const k2_y,
                          const double * const k2_s,
                          const double * const k3_x,
                          const double * const k3_y,
                          const double * const k3_s,
                          const double * const k4_x,
                          const double * const k4_y,
                          const double * const k4_s,
                          const double * const k5_x,
                          const double * const k5_y,
                          const double * const k5_s,
                          const double * const k6_x,
                          const double * const k6_y,
                          const double * const k6_s,
                          const double * const k7_x,
                          const double * const k7_y,
                          const double * const k7_s,
                          const double * const b)
{
  int index = threadIdx.x;
  out_x[index] = b[0]*k1_x[index] + b[1]*k2_x[index] + b[2]*k3_x[index]
                 + b[3]*k4_x[index] + b[4]*k5_x[index] + b[5]*k6_x[index]
                 + b[6]*k7_x[index];
  out_y[index] = b[0]*k1_y[index] + b[1]*k2_y[index] + b[2]*k3_y[index]
                 + b[3]*k4_y[index] + b[4]*k5_y[index] + b[5]*k6_y[index]
                 + b[6]*k7_y[index];
  out_s[0] = b[0]*k1_s[0] + b[1]*k2_s[0] + b[2]*k3_s[0] + b[3]*k4_s[0]
                 + b[4]*k5_s[0] + b[5]*k6_s[0] + b[6]*k7_s[0];
  out_s[1] = b[0]*k1_s[1] + b[1]*k2_s[1] + b[2]*k3_s[1] + b[3]*k4_s[1]
                 + b[4]*k5_s[1] + b[5]*k6_s[1] + b[6]*k7_s[1];
}

std::vector<triple>
dormand_prince( const double * const x,
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

  double *zero, *zer2;
  checkError(cudaMalloc(&zero,sizeof(double)*size));
  checkError(cudaMemset(zero,0,sizeof(double)*size));
  checkError(cudaMalloc(&zer2,sizeof(double)*2));
  checkError(cudaMemset(zer2,0,sizeof(double)*2));

  double *k1_x, *k1_y, *k1_s, *k2_x, *k2_y, *k2_s, *k3_x, *k3_y, *k3_s;
  double *k4_x, *k4_y, *k4_s, *k5_x, *k5_y, *k5_s, *k6_x, *k6_y, *k6_s;
  double *k7_x, *k7_y, *k7_s, *err_x, *err_y, *err_s;
  checkError(cudaMalloc(&k1_x,sizeof(double)*size));
  checkError(cudaMalloc(&k1_y,sizeof(double)*size));
  checkError(cudaMalloc(&k1_s,sizeof(double)*2));
  checkError(cudaMalloc(&k2_x,sizeof(double)*size));
  checkError(cudaMalloc(&k2_y,sizeof(double)*size));
  checkError(cudaMalloc(&k2_s,sizeof(double)*2));
  checkError(cudaMalloc(&k3_x,sizeof(double)*size));
  checkError(cudaMalloc(&k3_y,sizeof(double)*size));
  checkError(cudaMalloc(&k3_s,sizeof(double)*2));
  checkError(cudaMalloc(&k4_x,sizeof(double)*size));
  checkError(cudaMalloc(&k4_y,sizeof(double)*size));
  checkError(cudaMalloc(&k4_s,sizeof(double)*2));
  checkError(cudaMalloc(&k5_x,sizeof(double)*size));
  checkError(cudaMalloc(&k5_y,sizeof(double)*size));
  checkError(cudaMalloc(&k5_s,sizeof(double)*2));
  checkError(cudaMalloc(&k6_x,sizeof(double)*size));
  checkError(cudaMalloc(&k6_y,sizeof(double)*size));
  checkError(cudaMalloc(&k6_s,sizeof(double)*2));
  checkError(cudaMalloc(&k7_x,sizeof(double)*size));
  checkError(cudaMalloc(&k7_y,sizeof(double)*size));
  checkError(cudaMalloc(&k7_s,sizeof(double)*2));
  checkError(cudaMalloc(&err_x,sizeof(double)*size));
  checkError(cudaMalloc(&err_y,sizeof(double)*size));
  checkError(cudaMalloc(&err_s,sizeof(double)*2));

  double *a1, *a2, *a3, *a4, *a5, *a6, *a7;
  checkError(cudaMalloc(&a1,sizeof(double)*7));
  checkError(cudaMalloc(&a2,sizeof(double)*7));
  checkError(cudaMalloc(&a3,sizeof(double)*7));
  checkError(cudaMalloc(&a4,sizeof(double)*7));
  checkError(cudaMalloc(&a5,sizeof(double)*7));
  checkError(cudaMalloc(&a6,sizeof(double)*7));
  checkError(cudaMalloc(&a7,sizeof(double)*7));

  // Butcher Tableau constants

  a1[0] = 1/5; a1[1] = a1[2] = a1[3] = a1[4] = a1[5] = a1[6] = 0;
  a2[0] = 3/40; a2[1] = 9/40; a2[2] = a2[3] = a2[4] = a2[5] = a2[6] = 0;
  a3[0] = 44/45; a3[1] = -56/15; a3[2] = 32/9; a3[3] = a3[4] = a3[5] = a3[6] =0;
  a4[0] = 19372/6561; a4[1] = -25360/2187; a4[2] = 64448/6561;
    a4[3] = -212/729; a4[4] = a4[5] = a4[6] = 0;
  a5[0] = 9017/3168; a5[1] = -355/33; a5[2] = 46732/5247; a5[3] = 49/176;
    a5[4] = -5103/18656; a5[5] = a5[6] = 0;
  a6[0] = 35/384; a6[1] = 0; a6[2] = 500/1113; a6[3] = 125/192;
    a6[4] = -2187/6784; a6[5] = 11/84; a6[6] = 0;
  a7[0] = 5179/57600; a7[1] = 0; a7[2] = 7571/16695; a7[3] = 393/640;
    a7[4] = -92097/339200; a7[5] = 187/2100; a7[6] = 1/40;

  double t = t_start;

  int total_points = (int)((t_end - t_start)/h + 1);
  int sampling_rate = total_points / save;

  if ( sampling_rate <= 0 )
  {
    sampling_rate = 1;
  }
  int count = 0;

  cublasHandle_t handle;

  double error = 100;
  while(t <= t_end)
  {
    if(count % sampling_rate == 0)
    {
      SNAPSHOT(store,update_x,update_y,update_s,size);

      printf("current time: %f, final time: %f",t,t_end);
      std::cout << std::endl;
    }

    while(error > tolerance)
    {
      // k1 Step
      force(k1_x,k1_y,k1_s,update_x,update_y,update_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          zero,zero,zer2,zero,zero,zer2,zero,zero,zer2,zero,zero,zer2,
          zero,zero,zer2,a1,h);

      // k2 Step
      force(k2_x,k2_y,k2_s,temp_x,temp_y,temp_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          k2_x,k2_y,k2_s,zero,zero,zer2,zero,zero,zer2,zero,zero,zer2,
          zero,zero,zer2,a2,h);

      // k3 Step
      force(k3_x,k3_y,k3_s,temp_x,temp_y,temp_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          k2_x,k2_y,k2_s,k3_x,k3_y,k3_s,zero,zero,zer2,zero,zero,zer2,
          zero,zero,zer2,a3,h);

      // k4 Step
      force(k4_x,k4_y,k4_s,temp_x,temp_y,temp_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          k2_x,k2_y,k2_s,k3_x,k3_y,k3_s,k4_x,k4_y,k4_s,zero,zero,zer2,
          zero,zero,zer2,a4,h);

      // k5 Step
      force(k5_x,k5_y,k5_s,temp_x,temp_y,temp_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          k2_x,k2_y,k2_s,k3_x,k3_y,k3_s,k4_x,k4_y,k4_s,k5_x,k5_y,k5_s,
          zero,zero,zer2,a5,h);

      // k6 Step
      force(k6_x,k6_y,k6_s,temp_x,temp_y,temp_s,delta,p);
      construct_next_step<<<1,size>>>
          (temp_x,temp_y,temp_s,update_x,update_y,update_s,k1_x,k1_y,k1_s,
          k2_x,k2_y,k2_s,k3_x,k3_y,k3_s,k4_x,k4_y,k4_s,k5_x,k5_y,k5_s,
          k6_x,k6_y,k6_s,a6,h);

      // k7 Step
      force(k7_x,k7_y,k7_s,temp_x,temp_y,temp_s,delta,p);

      construct_error_step<<<1,size>>>
          (err_x,err_y,err_s,k1_x,k1_y,k1_s,k2_x,k2_y,k2_s,k3_x,k3_y,k3_s,
          k4_x,k4_y,k4_s,k5_x,k5_y,k5_s,k6_x,k6_y,k6_s,k7_x,k7_y,k7_s,a6);
      /*
      checkError(cudaMemcpy(heun,update_x,sizeof(double)*size,cudaD2D));
      checkError(cudaMemcpy(heun+size,update_y,sizeof(double)*size,cudaD2D));
      checkError(cudaMemcpy(heun+size*2,update_s,sizeof(double)*2,cudaD2D));

      // Error Handling
      vector_subtraction<<<1,size>>>(vec,euler,heun);
      cublasCreate(&handle);
      cublasSnrm2(handle,2*size+2,vec,1,&error); // Euclidean Norm

      if(error > tolerance)
        h /= 2; */
    }

    ++count;
    t+=h;
  }

  SNAPSHOT(store,update_x,update_y,update_s,size);

  // Wall of Destruction

  checkError(cudaFree(temp_x)); checkError(cudaFree(update_x));
  checkError(cudaFree(temp_y)); checkError(cudaFree(update_y));
  checkError(cudaFree(temp_s)); checkError(cudaFree(update_s));
  checkError(cudaFree(k1_x)); checkError(cudaFree(k5_x));
  checkError(cudaFree(k1_y)); checkError(cudaFree(k5_x));
  checkError(cudaFree(k1_s)); checkError(cudaFree(k5_x));
  checkError(cudaFree(k3_x)); checkError(cudaFree(k6_x));
  checkError(cudaFree(k3_y)); checkError(cudaFree(k6_x));
  checkError(cudaFree(k3_s)); checkError(cudaFree(k6_x));
  checkError(cudaFree(k2_x)); checkError(cudaFree(k7_x));
  checkError(cudaFree(k2_y)); checkError(cudaFree(k7_x));
  checkError(cudaFree(k2_s)); checkError(cudaFree(k7_x));
  checkError(cudaFree(k4_x));
  checkError(cudaFree(k4_y));
  checkError(cudaFree(k4_s));
  checkError(cudaFree(a1)); checkError(cudaFree(a4)); checkError(cudaFree(a6));
  checkError(cudaFree(a2)); checkError(cudaFree(a5)); checkError(cudaFree(a6));
  checkError(cudaFree(a3)); checkError(cudaFree(zero));
  checkError(cudaFree(zer2)); checkError(cudaFree(err_x));
  checkError(cudaFree(err_y)); checkError(cudaFree(err_s));

  return store;
}
