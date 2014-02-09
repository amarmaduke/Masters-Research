#include "utility.h"

namespace
{
  double * kahan_summation( int count, int size, double * constants,
                            double ** vectors)
  {
    double * sum = new double[count]();
    double * correction = new double[count]();

    for(int k = 0; k < count; ++k)
    {
      for(int i = 0; i < size; ++i)
      {
        double y = constants[k]*vectors[k][i] - correction[k];
        double t = sum[k] + y;
        correction[k] = (t - sum[k]) - y;
        sum[k] = t;
      }
    }

    delete[] correction;
    return sum;
  }
}

bool util::test::linc(int seed, bool verbose, int max_vector_count,
                      int max_vector_size, int max_int)
{
  srand(seed);
  int N = (rand() % max_vector_count)+1, s = (rand() % max_vector_size)+1;

  double **h_v, *c;
  h_v = new double*[N];
  c = new double[N];

  for(int i = 0; i < N; ++i)
  {
    h_v[i] = new double[s];
  }

  for(int i = 0; i < N; ++i)
  {
    c[i] = (rand() % max_int) + 1;
    for(int j = 0; j < s; ++j)
    {
      h_v[i][j] = rand() % max_int;
    }
  }

  thrust::host_vector<thrust::device_ptr<double> > vecs;
  thrust::host_vector<double> con;

  for(int i = 0; i < N; ++i)
  {
    double *d;
    cudaMalloc(&d,sizeof(double)*s);
    cudaMemcpy(d,h_v[i],sizeof(double)*s,cudaMemcpyHostToDevice);
    vecs.push_back(thrust::device_pointer_cast(d));

    con.push_back(c[i]);
  }

  double *h_result = kahan_summation(N, s, c, h_v);

  cublasHandle_t hand;
  cublasCreate(&hand);

  thrust::device_ptr<double> d_result = util::linc(hand,s,con,vecs);
  double *t_result;
  t_result = new double[s];
  cudaMemcpy(t_result,d_result.get(),sizeof(double)*s,cudaMemcpyDeviceToHost);

  double max_rel_error = 0, max_abs_error = 0;
  double avg_rel_error = 0, avg_abs_error = 0;

  for(int i = 0; i < s; ++i)
  {
    double abs_error = abs(h_result[i] - t_result[i]);
    double rel_error = abs_error / h_result[i];
    max_rel_error = max(max_rel_error,rel_error);
    max_abs_error = max(max_abs_error,abs_error);
    avg_rel_error = (rel_error + i*avg_rel_error)/(i+1);
    avg_abs_error = (abs_error + i*avg_abs_error)/(i+1);
  }
  if(verbose)
  {
    std::cout << "Max Relative Error: " << max_rel_error << std::endl;
    std::cout << "Max Absolute Error: " << max_abs_error << std::endl;
    std::cout << "Avg Relative Error: " << avg_rel_error << std::endl;
    std::cout << "Avg Absolulte Error: " << avg_abs_error << std::endl;
  }

  if(max_rel_error > .001)
    return false;
  else
    return true;
  //std::cout << "Test complete" << std::endl;
}
