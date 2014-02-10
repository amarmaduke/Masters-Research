#include <thrust/device_ptr.h>
#include <thrust/transform.h>

#include "parameter.h"
#include "../utility/utility.h"

#ifndef FORCE_H
#define FORCE_H

struct force_functor
{
  parameter state;

  force_functor(parameter p) : state(p) { };

  thrust::device_ptr<double> operator() (double t,thrust::device_ptr<double> y);
};

#endif
