#include "dorpi.h"
#include "../simulation/force.h"

#define cudaH2D cudaMemcpyHostToDevice
#define cudaD2H cudaMemcpyDeviceToHost
#define cudaD2D cudaMemcpyDeviceToDevice

void SNAPSHOT(std::vector<triple>& save, double *x, double *y, double *s, int size)
{
  double *first = new double[size];
  double *second = new double[size];
  double *third = new double[2];

  check_error(cudaMemcpy(first,x,sizeof(double)*size,cudaD2H));
  check_error(cudaMemcpy(second,y,sizeof(double)*size,cudaD2H));
  check_error(cudaMemcpy(third,s,sizeof(double)*2,cudaD2H));

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

double * linc_s(cublasHandle_t handle, int count, int size,
              std::vector<double>& constants, std::vector<double *>& vectors)
{
  // Preconditions
  assert(constants.size() == vectors.size());
  count = count < 1 || count > vectors.size() ? vectors.size() : count;

  std::vector<int> iter(count);
  std::vector<double *> vectors_copy(vectors.size());
  std::vector<cudaStream_t> streams(count);

  double *buffer;
  check_error(cudaMalloc(&buffer, sizeof(double)*size*count));
  for(int i = 0; i < count; ++i)
  {
    // Setup iterator
    iter.push_back(i);

    // Setup Streams
    cudaStream_t stream;
    check_error(cudaStreamCreate(&stream));
    streams.push_back(stream);

    // Grab vector
    double *temp = buffer + i*size;
    cublasSetStream(handle,stream);
    cublasDcopy(handle, size, vectors[i], 1, temp, 1);
    vectors_copy.push_back(temp);
  }
  cudaDeviceSynchronize();

  // Compute linear compination
  int N = count, m, j = 0;
  while(N > 1)
  {
    m = N % 2 == 0 ? N : N - 1;
    for(int i = 0; i < m; i+=2, ++j)
    {
      double alpha = constants[iter[i+1]]/constants[iter[i]];
      cublasSetStream(handle,streams[i]);
      cublasDaxpy(handle, count, &alpha, vectors_copy[iter[i+1]], 1,
                  vectors_copy[iter[i]], 1);
      iter[j] = iter[i];
    }
    cudaDeviceSynchronize();
    if(N % 2 != 0)
      iter[j] = iter[m];
    N = (N+1)/2;
    j = 0;
  }

  double a = constants[iter[0]];
  cublasDscal(handle, count, &a, vectors[0], 1);

  // Grab output and clean
  double * out = vectors[0];
  check_error(cudaFree(buffer));
  for(int i = 0; i < count; ++i)
  {
    cudaStreamDestroy(streams[i]);
  }

  return out;
}
