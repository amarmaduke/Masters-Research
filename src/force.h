#ifndef _FORCE_H_
#define _FORCE_H_

#include "parameter.h"

void force(	double * const out_x,
						double * const out_y,
						double * const out_s,
						const double * const in_x,
						const double * const in_y,
						const double * const in_s,
						const parameter h_p,
						const parameter * const d_p);

#endif
