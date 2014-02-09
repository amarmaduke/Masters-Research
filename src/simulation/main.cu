#include <iostream>
#include <vector>
#include "euler.h"
#include "parameter.h"
#include "dorpi.h"
#include <stdio.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"

#define cudaH2D cudaMemcpyHostToDevice
#define cudaD2H cudaMemcpyDeviceToHost
#define cudaD2D cudaMemcpyDeviceToDevice

using namespace std;

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

void print(vector<triple> in, int s)
{
  for(int i = 0; i < in.size(); ++i)
  {
    triple p = in[i];

    for(int j = 0; j < s; ++j)
    {
      printf("(%.10f,%.10f) ",p.first[j],p.second[j]);
      cout << endl;
    }
    printf("(%.10f,%.10f)",p.third[0],p.third[1]);
    cout << endl << endl;
  }
  cout << endl << endl;
}

int main()
{
	test();
	/*
  double t = 0;
  double step = .0001;
  double t_end = .05;
  parameter p;
  int size = p.m*p.n;

  double* initial_x = new double[size];
  double* initial_y = new double[size];
  double* initial_s = new double[2];
  double* initial_delta = new double[p.m];
  for(int j = 0; j < p.m; ++j)
  {
    initial_delta[j] = j;
    for(int i = 0; i < p.n; ++i)
    {
      int index = i + p.n*j;
      initial_x[index] = initial_delta[j];
      initial_y[index] = i+1;
    }
  }
  initial_s[0] = -5;
  initial_s[1] = 12;

  double *x, *y, *s, *delta;

  checkError(cudaMalloc(&x,sizeof(double)*size));
  checkError(cudaMalloc(&y,sizeof(double)*size));
  checkError(cudaMalloc(&s,sizeof(double)*2));
  checkError(cudaMalloc(&delta,sizeof(double)*p.m));

  checkError(cudaMemcpy(x,initial_x,sizeof(double)*size,cudaH2D));
  checkError(cudaMemcpy(y,initial_y,sizeof(double)*size,cudaH2D));
  checkError(cudaMemcpy(s,initial_s,sizeof(double)*2,cudaH2D));
  checkError(cudaMemcpy(delta,initial_delta,sizeof(double)*p.m,cudaH2D));

  time_t start = time(0);
  vector<triple> grid
    = eulers_method(x,y,s,delta,t,t_end,step,2,p);
  time_t end = time(0);
  double time = difftime(end,start);

  print(grid,size);
  cout << "time: " << time << endl;

  checkError(cudaFree(x));
  checkError(cudaFree(y));
  checkError(cudaFree(s));
  checkError(cudaFree(delta));
  return 0;
	*/
}

