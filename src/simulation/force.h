#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

#include "parameter.h"
#include "../utility/utility.h"

#ifndef FORCE_H
#define FORCE_H

struct force_functor
{
  parameter state;

  force_functor(parameter p) : state(p) { };

  void operator() (const thrust::device_vector< double > &x,
                   thrust::device_vector< double > &dxdt,
                   const double dt);
};

#endif
