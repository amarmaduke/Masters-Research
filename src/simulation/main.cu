#include <iostream>
#include <vector>
#include <stdio.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "force.h"
#include "../solver/dorpi.h"

using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;

struct test_functor
{

  void operator() (const state_type &x, state_type &dxdt, const value_type dt)
  {
    thrust::copy(x.begin(), x.end(), dxdt.begin());
  }

};

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
		//std::cout << "t: " << t << std::endl;
	}
};

void ode_test()
{
  parameter p;
  int size = p.m*p.n;
  thrust::host_vector< value_type > init(2*size+2);
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

  force_functor F(p);
  test_observer O(g_time_point,g_save);

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
	{ // Burned run
  force_functor F(p);
	std::vector<double> v;
	std::vector<double*> vp;

  test_observer O(v,vp);

  state_type x = init;

	integrate_const(make_controlled(1.0e-11, 1.0e-11, stepper_type()),
                  F, x, 0.0, 10.0, 1., O);
	}

	cudaEvent_t start, stop;
	float time;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0);
	for(int i = 0; i < 0; ++i)
	{
		force_functor F(p);
		std::vector<double> v;
		std::vector<double*> vp;

  	test_observer O(v,vp);

  	state_type x = init;
		std::cout << "i: " << i << std::endl;
		integrate_const(make_controlled(1.0e-11, 1.0e-11, stepper_type()),
                  F, x, 0.0, 10.0, 1., O);
	}
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);


	/*
	for(int i = 0; i < g_save.size(); ++i)
	{
		double* s = g_save[i];
		std::cout << "t: " << g_time_point[i] << std::endl;
		for(int j = 0; j < 2*size+2; ++j)
		{
			std::cout << s[j] << " ";
		}
		std::cout << std::endl;
	}*/

	std::cout << "Total time: " << time << " ms, " << time/1000. << " s" << std::endl;
}

int main()
{
  ode_test();

  /*
  double t = 0;
  double step = .01;
  double t_end = .05;
  parameter p;
  int size = p.m*p.n;

  thrust::device_ptr<double> init = thrust::device_malloc<double>(2*size+2);
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

  force_functor F(p);
  dorpi::options o;
  o.t_start = t; o.t_end = t_end; o.requested_h = step;
  o.tolerance = .00001;
	o.save_count = 10;

	util::print(init,20);
	util::print(d_delta,p.m);

  time_t start = time(0);
  thrust::host_vector<thrust::device_ptr<double> > grid
      = dorpi::solve(F, t, init, 2*size+2, o);
  time_t end = time(0);
  double time = difftime(end,start);

	std::cout << "Grid:" << std::endl;
  for(int i = 0; i < grid.size(); ++i)
  {
		thrust::device_ptr<double> d_r = grid[i];
		util::print(d_r,2*size+2);
	}
  std::cout << "time: " << time << std::endl;

  thrust::device_free(init);
  thrust::device_free(d_delta);
  */
  return 0;
}
