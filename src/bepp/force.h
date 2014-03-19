#include <limits>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

#include "parameter.h"
#include "defs.h"

#ifndef FORCE_H
#define FORCE_H

struct force_functor
{
  parameter& state;

  force_functor(parameter& p) : state(p) { };

  void operator() (const vector_type &x,
                   vector_type &dxdt,
                   const value_type dt);
};

#endif
