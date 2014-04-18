#ifndef DEFS_H
#define DEFS_H

#include <sundials/sundials_types.h>

#define ZERO RCONST(0)
#define ONE RCONST(1)
#define TWO RCONST(2)
#define FOUR RCONST(4)
#define EIGHTH RCONST(.125)
#define PI RCONST(3.141592653589793238463)

// Thread block size for n-body kernel
// Pick multiple of 32, less than or equal to 1024
#define K 32

// Number of simulations to compute in parallel for pulloff profile
#define SIM_COUNT 10

// Utilities
template<typename T, typename U>
inline T* as(U* u) { return dynamic_cast<T*>(u); }

template<typename T, typename U>
inline const T* as(const U* u) { return dynamic_cast<const T*>(u); }

template<typename T, typename U>
inline T& as(U& u) { return dynamic_cast<T&>(u); }

template<typename T, typename U>
inline const T& as(const U& u) { return dynamic_cast<const T&>(u); }

#endif
