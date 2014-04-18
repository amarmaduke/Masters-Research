
#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cmath>

#include <nvector/nvector_serial.h>  /* serial N_Vector types, */
#include <sundials/sundials_types.h> /* definition of realtype */

#include "parameter.h"
#include "defs.h"

extern "C" int force(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd_fiber,
                    std::vector<int>& nhbd_partner,
                    std::vector<bool>& mask,
                    parameter& param);

struct force_wrapper
{
  parameter& params;
  /*
    Two parallel vectors contain the neighborhood list for van der Waals. The
    first vector, nhbd_fiber will always be an index to a particle on a fiber.
    The second vector, nhbd_partner can either be a particle on a fiber, the
    upper substrate, or the lower substrate. The bitset vector, mask, will
    determine which indexing scheme the nhbd_partner index belongs to. Every
    two boolean values will be used to determine the correct index scheme.

    True, True (or 1, 1) will be a fiber particle,
    True, False (or 1, 0) will be an upper substrate particle,
    False, True (or 0, 1) will be a lower substrate particle,
    False, False (or 0, 0) is undefined and will cause undefined behavior.
  */
  std::vector<int> nhbd_fiber;
  std::vector<int> nhbd_partner;
  std::vector<bool> mask;

  N_Vector& u_ref;

  force_wrapper(parameter& p, N_Vector& ref, N_Vector& init)
  : params(p), u_ref(ref), nhbd_fiber(), nhbd_partner(), mask()
  {
    generate_nhbd(init, nhbd_fiber, nhbd_partner, mask, params);
  }
};

#endif
