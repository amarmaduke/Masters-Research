#include <thrust/device_vector.h>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

#ifndef DEFS_H
#define DEFS_H

// double or float
typedef
  double
  value_type;

// device_vetor, host_vector, std::vector, ...
typedef
  thrust::device_vector<value_type>
  vector_type;

// Any explicit solver in boost::...::odeint
typedef
  runge_kutta_dopri5< vector_type, value_type, vector_type, value_type >
  stepper_type;

// Thread block size for n-body kernel
#define K 32

const double PI = 3.141592653589793238463;

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
