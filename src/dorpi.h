#ifndef _DORMAND_PRINCE_H_
#define _DORMAND_PRINCE_H_

#include "parameter.h"
#include <vector>
#include <utility>

std::vector<triple>
dormand_prince( const double* const x,
                const double* const y,
                const double* const s,
                const double* const delta,
                double t_start,
                double t_end,
                double h,
                int save,
                double tolerance,
                parameter p);

#endif
