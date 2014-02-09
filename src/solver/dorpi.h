#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"

#include "../utility/utility.h"
#include "../simulation/parameter.h"

#ifndef _DORMAND_PRINCE_H_
#define _DORMAND_PRINCE_H_

namespace dorpi
{

struct options
{
  double t_start;
  double t_end;
  double requested_h;
  double tolerance;
  int save_count;
};

template<typename F>
thrust::host_vector<thrust::device_ptr<double> >
solve(F& functor, thrust::device_ptr<double>& t_init,
      thrust::device_ptr<double>& y_init, int ptr_size, const options& o);

#include "dorpi.imp"

} // dorpi

#endif
