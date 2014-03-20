#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

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
		if(p.have_delta)
			d_delta[j] = p.delta[j];
		else
    	d_delta[j] = j;
    for(int i = 0; i < p.n; ++i)
    {
      int index = i + p.n*j;
			if(p.have_init)
			{
				init[index] = p.init[index];
				init[index+size] = p.init[index+size];
			}else
			{
      	init[index] = d_delta[j];
      	init[index+size] = i+1;
    	}
		}
  }
	if(p.have_init)
	{
		init[0+2*size] = p.init[0+2*size];
		init[1+2*size] = p.init[1+2*size];
		free(p.delta);
	}
  else
	{
		init[0+2*size] = p.sub_x;
  	init[1+2*size] = p.sub_y;
  }
	p.delta = d_delta.get();

  cudaEvent_t start, stop;
  float timer;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);

  std::vector< value_type > v;
  std::vector< value_type* > vp;
  std::vector< value_type > times;
  times.push_back(0);
  times.push_back(9);
  times.push_back(10);

	force_functor F(p);
  observer O(v,vp);

	value_type dt = .1;
  value_type t0 = 0;
	//integrate_times(make_controlled(p.abstol, p.reltol, stepper_type()),
  //              F, init, times.begin(), times.end(), dt, O);
	integrate_n_steps(make_controlled(p.abstol, p.reltol, stepper_type()),
								F, init, t0, dt, 10000, O);

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  for(int i = 0; i < vp.size(); ++i)
  {
		json::Array* temp = new json::Array;
    value_type* s = vp[i];
		std::stringstream ss;
		ss << "t:" << v[i];
		std::string temp2 = ss.str();
    for(int j = 0; j < 2*size+2; ++j)
    {
			json::Number* num = new json::Number;
			num->val = s[j];
			temp->push_back(num);
    }
		obj[temp2] = temp;
  }

	std::cout << "time: " << timer << " ms " << (timer/1000) << " s" << std::endl;

	std::ofstream File("output.json");
	json::print(File,obj);
	std::ofstream File2("/home/marmaduke/output.json");
	json::print(File2,obj);
}

int main()
{
  json::Object obj = json::parse(std::cin,10);
  ode_test(obj);
  return 0;
}
