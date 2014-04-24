
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
  realtype LJ_f2f_c, LJ_f2l_c, LJ_f2u_c;

  force_wrapper(parameter& p, N_Vector& init)
  : params(p), nhbd_fiber(), nhbd_partner(), mask()
  {
    generate_nhbd(init, nhbd_fiber, nhbd_partner, mask, params);

    realtype p1, p2, p4, p8, p7, p13;
    p1 = params.sigma/params.rcut;
    p2 = p1*p1; p4 = p2*p2; p8 = p4*p4;
    p7 = p4*p2*p1; p13 = p8*p4*p1;
    LJ_f2u_c = -(RCONST(12)*params.epsilon_top/params.sigma)*(p13 - p7);
    LJ_f2l_c = -(RCONST(12)*params.epsilon_bottom/params.sigma)*(p13 - p7);
    LJ_f2f_c = -(RCONST(12)*params.epsilon/params.sigma)*(p13 - p7);
  }
};

#endif
