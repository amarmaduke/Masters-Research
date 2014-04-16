
#ifndef FORCE_H
#define FORCE_H

#include <iostream>

#include <nvector/nvector_serial.h>  /* serial N_Vector types, */
#include <sundials/sundials_types.h> /* definition of realtype */

#include "parameter.h"

struct force_wrapper
{
  parameter& params;
  int* nhbd_i;
  int* nhbd_j;

  force_wrapper(parameter& p) : params(p) { }
};

extern "C" int force(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#endif
