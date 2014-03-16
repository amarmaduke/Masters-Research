#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "force.h"
#include "../solver/dorpi.h"
#include "../json/json.h"
#include "main.h"

using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::host_vector< value_type > state_type;

std::vector<double> g_time_point;
std::vector<double*> g_save;

struct test_observer
{
	std::vector<double>& time_point;
	std::vector<double*>& save;

	test_observer(std::vector<double>& t,
			std::vector<double* >& s) : time_point(t), save(s) { }

  template<typename State >
  void operator() (const State &x, value_type t )
  {
		time_point.push_back(t);
		double* s = new double[x.size()];
		thrust::copy(x.begin(),x.end(),s);
  	save.push_back(s);
		std::cout << "t: " << t << std::endl;
	}
};

struct force_test
{
	void operator() (const thrust::host_vector<double>& x,
										thrust::host_vector<double>& dxdt,
										const double t)
	{
		for(int i = 0; i < x.size(); ++i)
		{
			if(x[i] > 0)
			{
				for(int j = 0; j < 2; ++j)
				{
					dxdt[i] = 1/x[i];
				}
			}else
			{
				dxdt[i] = -x[i];
			}
		}
	}
};

void ode_test(json::Object& obj)
{
	cudaStream_t s;
	cudaStreamCreate(&s);

  parameter p(obj);
  int size = p.m*p.n;
	value_type * temp_init;
	cudaHostAlloc(&temp_init,sizeof(value_type)*2+2,cudaHostAllocPortable);
	thrust::host_vector<value_type> init(temp_init,temp_init+size*2+2);
  thrust::device_ptr<double> d_delta = thrust::device_malloc<double>(p.m);
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
  init[0+2*size] = -5;
  init[1+2*size] = 12;
  p.delta = d_delta.get();

  typedef
    runge_kutta_dopri5< state_type, value_type, state_type, value_type >
    stepper_type;

  state_type x = init;

	std::vector<double> time_range;
	//state_type temp(2*size + 2);

/*
	temp.begin();
	F(x,temp,0);
	double t_d = temp[2*size];
	printf("%.12f\n",t_d);
	for(int i = 0; i < size; ++i)
	{
		t_d = temp[i];
		printf("%4.12f\n",t_d);
	}
	t_d = temp[2*size+1];
	printf("%4.12f\n",t_d);
	for(int i = size; i < 2*size; ++i)
	{
		t_d = temp[i];
		printf("%4.12f\n",t_d);
	}
  */
	/*{ // Burned run
  force_functor F(p);
	std::vector<double> v;
	std::vector<double*> vp;

  test_observer O(v,vp);

  state_type x = init;

	integrate_const(make_controlled(1.0e-11, 1.0e-11, stepper_type()),
                  F, x, 0.0, 10.0, 1., O);
	}*/

	cudaEvent_t start, stop;
	float time;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0);
	for(int i = 0; i < 1; ++i)
	{
		force_functor F(p,s);
		std::vector<double> v;
		std::vector<double*> vp;

		test_observer O(v,vp);

  	state_type x = init;
		std::cout << "i: " << i << std::endl;
		integrate_const(make_controlled(1.0e-11, 1.0e-11, stepper_type()),
                  F, x, 0.0, 100.0, 20., O);
		
		for(int i = 0; i < vp.size(); ++i)
		{
			double* s = vp[i];
			std::cout << "t: " << v[i] << std::endl;
			for(int j = 0; j < 2*size+2; ++j)
			{
				printf("%3.11f ",s[j]);
			}
			std::cout << std::endl;
		}
	}
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	cudaStreamDestroy(s);

	std::cout << "Total time: " << time << " ms, " << time/1000. << " s" << std::endl;
}

