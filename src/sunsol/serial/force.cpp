
#include "force.h"

extern "C" int force(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  force_wrapper wp = *((force_wrapper*)(user_data));
  parameter p = wp.params;
  std::cout << p.n << " " << p.m << std::endl;
  N_VAddConst(y, 0, ydot);

  return 0;
}
