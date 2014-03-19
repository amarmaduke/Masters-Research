#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <thrust/device_vector.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "force.h"
#include "defs.h"
#include "../json/json.h"

using namespace boost::numeric::odeint;

struct observer
{
  std::vector<value_type>& time_point;
  std::vector<value_type*>& save;

  observer(std::vector<value_type>& t,
      std::vector<value_type* >& s) : time_point(t), save(s) { }

  template<typename State >
  void operator() (const State &x, value_type t )
  {
    time_point.push_back(t);
    value_type* s = new value_type[x.size()];
    thrust::copy(x.begin(),x.end(),s);
    save.push_back(s);
    std::cout << "t: " << t << std::endl;
  }
};

void ode_test(json::Object& obj)
{
  parameter p(obj);
  int size = p.m*p.n;

  vector_type init(size*2+2);
  thrust::device_ptr< value_type > d_delta =
          thrust::device_malloc< value_type >(p.m);
  for(int j = 0; j < p.m; ++j)
  {
    d_delta[j] = j;
    for(int i = 0; i < p.n; ++i)
    {
      int index = i + p.n*j;
      init[index] = d_delta[j];
      init[index+size] = i+1;
    }
  }
  init[0+2*size] = p.sub_x;
  init[1+2*size] = p.sub_y;
  p.delta = d_delta.get();

  cudaEvent_t start, stop;
  float timer;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);

  force_functor F(p);
  std::vector< value_type > v;
  std::vector< value_type* > vp;
  std::vector< value_type > times;
  times.push_back(0);
  times.push_back(9);
  times.push_back(10);


  observer O(v,vp);

  std::cout << p.abstol << " " << p.reltol << std::endl;

  integrate_times(make_controlled(p.abstol, p.reltol, stepper_type()),
                F, init, times.begin(), times.end(), .01, O);

  for(int i = 0; i < vp.size(); ++i)
  {
    value_type* s = vp[i];
    std::cout << "t: " << v[i] << std::endl;
    for(int j = 0; j < 2*size+2; ++j)
    {
      printf("%3.11f ",s[j]);
    }
    std::cout << std::endl;
  }

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  std::cout << "Total time: " << timer << " ms, " << timer/1000. << " s" << std::endl;
}

int main()
{
  json::Object obj = json::parse(std::cin,10);
  ode_test(obj);
  return 0;
}
