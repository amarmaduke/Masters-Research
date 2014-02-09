#ifndef _EULER_H_
#define _EULER_H_

#include "parameter.h"
#include <vector>
#include <utility>

std::vector<triple>
eulers_method(const double* const x,
              const double* const y,
              const double* const s,
              const double* const delta,
              double t_start,
              double t_end,
              double h,
              int save,
              parameter p);

#endif
