#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"

#include <cstdio>
#include <cassert>
#include <typeinfo>

#include "def.h"

#ifndef _UTILITY_H_
#define _UTILITY_H_

#ifdef _ERROR_
  #define check_error(x) do { util::error_code e(x); \
                              util::gpu_assert(e); } while(0)
#else
  #define check_error(x) x
#endif

/*
  Collection of utility and helper functions.
*/

namespace util
{

struct error_code
{
  enum { CUDA, CUBLAS } tag;
  union {
    cudaError_t cuda;
    cublasStatus_t cublas;
  } type;

  error_code(cudaError_t e);
  error_code(cublasStatus_t e);

  operator const char*() const;
  operator int() const;
  bool operator!() const;
};

void gpu_assert(error_code);

thrust::device_ptr<double> linc(cublasHandle_t handle, int size,
                    thrust::host_vector<double>& constants,
                    thrust::host_vector<thrust::device_ptr<double> >& vectors);


template<typename T>
struct pinned_ptr {

  __host__
  pinned_ptr() : data(0) {};

  __host__
  pinned_ptr(T* host_ptr) : data(host_ptr) {};

  __host__
  T* get() const;

/*
  __host__
  pinned_ptr& operator=(const pinned_ptr& ptr);
*/

  private:
   T* data;
};

template<typename T>
__host__
pinned_ptr<T> inline pinned_malloc(size_t);

template<typename T>
__host__ inline
void pinned_free(pinned_ptr<T>);

template<typename T>
__host__ inline
void deep_copy(thrust::device_ptr<T> dest, thrust::device_ptr<T> src, size_t N);

template<typename T>
__host__ inline
void deep_copy(thrust::device_ptr<T> dest, T* src, size_t);

template<typename T>
__host__ inline
void deep_copy(T* dest, thrust::device_ptr<T> src, size_t);

template<typename T>
__host__ inline
void deep_copy( thrust::device_ptr<T> dest, thrust::device_ptr<T> src, size_t,
                cudaStream_t);

template<typename T>
__host__ inline
void deep_copy( thrust::device_ptr<T> dest, pinned_ptr<T> src, size_t,
                cudaStream_t);

template<typename T>
__host__ inline
void deep_copy( pinned_ptr<T> dest, thrust::device_ptr<T> src, size_t,
                cudaStream_t);

namespace test {

bool linc(int seed, bool verbose, int max_vector_count, int max_vector_size,
          int max_int);

} // util::test

#include "utility.imp"

} // util

#endif //_UTILITY_H_
