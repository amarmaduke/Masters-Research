#include "utility.h"

static const char * cublasGetErrorString(cublasStatus_t error)
{
  switch (error)
  {
      case CUBLAS_STATUS_SUCCESS:
          return "CUBLAS_STATUS_SUCCESS";

      case CUBLAS_STATUS_NOT_INITIALIZED:
          return "CUBLAS_STATUS_NOT_INITIALIZED";

      case CUBLAS_STATUS_ALLOC_FAILED:
          return "CUBLAS_STATUS_ALLOC_FAILED";

      case CUBLAS_STATUS_INVALID_VALUE:
          return "CUBLAS_STATUS_INVALID_VALUE";

      case CUBLAS_STATUS_ARCH_MISMATCH:
          return "CUBLAS_STATUS_ARCH_MISMATCH";

      case CUBLAS_STATUS_MAPPING_ERROR:
          return "CUBLAS_STATUS_MAPPING_ERROR";

      case CUBLAS_STATUS_EXECUTION_FAILED:
          return "CUBLAS_STATUS_EXECUTION_FAILED";

      case CUBLAS_STATUS_INTERNAL_ERROR:
          return "CUBLAS_STATUS_INTERNAL_ERROR";
  }

  return "<unknown>";
}

util::error_code::error_code(cudaError_t e)
{
  this->tag = CUDA;
  this->type.cuda = e;
}

util::error_code::error_code(cublasStatus_t e)
{
  this->tag = CUBLAS;
  this->type.cublas = e;
}

util::error_code::operator const char*() const
{
  switch(this->tag)
  {
    case CUDA: return cudaGetErrorString(this->type.cuda);
    case CUBLAS: return cublasGetErrorString(this->type.cublas);
    default: return "Unknown Error Code";
  }
}

util::error_code::operator int() const
{
  switch(this->tag)
  {
    case CUDA: return this->type.cuda;
    case CUBLAS: return this->type.cublas;
    default: return -1;
  }
}

bool util::error_code::operator!() const{
  switch(this->tag)
  {
    case CUDA: return this->type.cuda == cudaSuccess ? true : false;
    case CUBLAS: return this->type.cublas == CUBLAS_STATUS_SUCCESS ?
          true : false;
    default: return true;
  }
}

void util::gpu_assert(error_code code)
{
  if (!code)
  {
    fprintf(stderr,"cuda_assert: %s %s %d\n",
      (const char*)code, __FILE__, __LINE__);
    exit(code);
  }
}

thrust::device_ptr<double> util::linc(cublasHandle_t handle, int size,
                    thrust::host_vector<double>& constants,
                    thrust::host_vector<thrust::device_ptr<double> >& vectors)
{
  // Preconditions
  assert(constants.size() == vectors.size());

  // Setup
  int count = constants.size();
  std::vector<int> iter(count);
  thrust::host_vector<thrust::device_ptr<double> > vectors_copy(vectors.size());

  thrust::device_ptr<double> buffer = thrust::device_malloc<double>(size*count);
  for(int i = 0; i < count; ++i)
  {
    // Setup iterator
    iter[i] = i;

    // Grab vector pointer
    vectors_copy[i] = buffer + size*i;
    util::deep_copy(vectors_copy[i],vectors[i],size);
  }

  // Compute linear combination
  int N = count, m, j = 0;
  while(N > 1)
  {
    m = N % 2 == 0 ? N : N - 1;
    for(int i = 0; i < m; i+=2, ++j)
    {
      double alpha = constants[iter[i+1]]/constants[iter[i]];

      check_error(cublasDaxpy(handle, size, &alpha,
                              vectors_copy[iter[i+1]].get(), 1,
                              vectors_copy[iter[i]].get(), 1));

      iter[j] = iter[i];
    }
    if(N % 2 != 0)
      iter[j] = iter[m];
    N = (N+1)/2;
    j = 0;
  }

  double a = constants[iter[0]];
  check_error(cublasDscal(handle, size, &a, vectors_copy[iter[0]].get(), 1));

  // Grab output and clean
  thrust::device_ptr<double> out = thrust::device_malloc<double>(size);
  util::deep_copy(out, vectors_copy[iter[0]], size);
  thrust::device_free(buffer);

  return out;
}
