#ifndef _DEF_H_
#define _DEF_H_

#define _ERROR_
#define _DEBUG_
#define cudaH2D cudaMemcpyHostToDevice
#define cudaD2H cudaMemcpyDeviceToHost
#define cudaD2D cudaMemcpyDeviceToDevice

const double PI = 3.141592653589793238463;

template<typename T, typename U>
inline T* as(U* u) { return dynamic_cast<T*>(u); }

template<typename T, typename U>
inline const T* as(const U* u) { return dynamic_cast<const T*>(u); }

template<typename T, typename U>
inline T& as(U& u) { return dynamic_cast<T&>(u); }

template<typename T, typename U>
inline const T& as(const U& u) { return dynamic_cast<const T&>(u); }

#endif // _DEF_H_
