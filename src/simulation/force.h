#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>

#include <cmath>

#include "parameter.h"
#include "../utility/def.h"
#include "../utility/utility.h"

#ifndef FORCE_H
#define FORCE_H

struct force_functor
{
  parameter& state;
	cudaStream_t& stream;

  force_functor(parameter& p, cudaStream_t& s) : state(p), stream(s) { };

  void operator() (const thrust::host_vector< double > &x,
                   thrust::host_vector< double > &dxdt,
                   const double dt);
};

#endif
