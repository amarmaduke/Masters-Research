#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "force.h"
#include "defs.h"
#include "../json/json.h"

using namespace boost::numeric::odeint;

value_type norm(vector_type& a, vector_type& b)
{
	thrust::plus<value_type> add;
	thrust::minus<value_type> sub;
	thrust::multiplies<value_type> mul;
	thrust::transform(a.begin(),a.end(),b.begin(),b.begin(),sub);
	thrust::transform(b.begin(),b.end(),b.begin(),b.begin(),mul);
	value_type result = thrust::reduce(b.begin(),b.end(),0.0,add);
	return sqrt(result);
}

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

struct pp_observer
{
  vector_type& prev;
  vector_type& curr;
	json::Object& obj;
	bool& swtch;
	int& count;

	pp_observer(vector_type& p, vector_type& c, json::Object& o, bool& s, int& co) 
		: prev(p), curr(c), obj(o), swtch(s), count(co) { }

	template<typename State >
	void operator() (const State&x, value_type t)
	{
		if(swtch)
		{
			thrust::copy(x.begin(),x.end(),prev.begin());
			swtch = !swtch;
		}else
		{
			thrust::copy(x.begin(),x.end(),curr.begin());
			swtch = !swtch;
			if(count != -1)
			{
				json::Array* a = new json::Array();
				for(int i = 0; i < curr.size(); ++i)
				{
					double x = curr[i];
					a->push_back(new json::Number(x));
				}
				std::stringstream ss;
				ss << "tq" << count;
				obj[ss.str()] = a;
				++count;
			}
		}
	}
};

void pulloff_profile(json::Object& obj)
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

  std::vector<value_type> times;
  vector_type prev(2*size+2);
	vector_type curr(2*size+2);
	times.push_back(9); times.push_back(10);

	bool swtch = false;
	int obs_count = -1;
	force_functor F(p);
	pp_observer O(prev,curr,obj,swtch,obs_count);
	vector_type V = vector_type(init);

  json::Array* grid = new json::Array();

  double linear_step = 5;
  int output_i = 1;
  double magnitude = 10;
  double lower_bound = 0;
  double upper_bound = 100;
  bool have_upper = false, have_lower = false;

  double dt = .1;
  double theta_begin = 15, theta_end = 45, theta_h = 5;

  for(int t = theta_begin; t <= theta_end; t+=theta_h)
  {
    bool found = false;
    have_lower = have_upper = false;

    if(t != theta_begin)
    {
      linear_step = 1;
    }
		std::cout << "Searching... theta = " << t << std::endl;

    while(not found)
    {
      bool done = false;
      int outcome = -1;
			
      double l = -magnitude*sin(PI*t/180);
      double m = magnitude*cos(PI*t/180);
			F.state.lambda = l;
			F.state.mu = m;
			while(not done)
      {
        integrate_times(make_controlled(p.abstol, p.reltol, stepper_type()),
          F, V, times.begin(), times.end(), dt, O);

				vector_type diff(2*size+2);
				thrust::minus<value_type> op;
				thrust::transform(curr.begin(),curr.end(),prev.begin(),diff.begin(),op);
				value_type equil = sqrt(thrust::inner_product(
								diff.begin(),diff.begin()+2*size,diff.begin(),0.));
    		value_type adhesion = sqrt(thrust::inner_product(
								diff.begin()+2*size,diff.end(),diff.begin()+2*size,0.));
    		value_type d = sqrt(l*l + m*m);
				std::cout << "equil: " << equil << " adhesion: " << adhesion << " d: " << d << std::endl;

        if(abs(adhesion - d) < p.abstol)
        { // Pulled off
          done = true;
          json::Array* temp = new json::Array();
          temp->push_back(new json::Number(t));
					temp->push_back(new json::Number(l)); 
					temp->push_back(new json::Number(m)); 
					temp->push_back(new json::Number(1));
          grid->push_back(temp);
          outcome = 1;
					std::cout << "Pulled off. l: " << l << " m: " << m << std::endl;
        }else if(equil < p.abstol and adhesion < p.abstol)
        { // Adhered
          done = true;
          json::Array* temp = new json::Array();
					temp->push_back(new json::Number(t));
          temp->push_back(new json::Number(l)); 
					temp->push_back(new json::Number(m)); 
					temp->push_back(new json::Number(0));
          grid->push_back(temp);
          outcome = 0;
					std::cout << "Adhered. l: " << l << " m: " << m << std::endl;
        }
      }

      if(have_upper and have_lower and abs(upper_bound - lower_bound) < .5)
      {
        found = true;
				V = vector_type(init);
      }

      if(not have_lower or not have_upper)
      {
        if(outcome == 0) // Adhered
        {
          lower_bound = magnitude;
          have_lower = true;
          if(not have_upper)
            magnitude = magnitude + linear_step;
        }else if(outcome == 1) // Pulled Off
        {
          upper_bound = magnitude;
          have_upper = true;
          if(not have_lower)
            magnitude = magnitude - linear_step;
        }
      }

      if(have_lower and have_upper)
      {
        if(outcome == 1) //Pulled off
          lower_bound = magnitude;
        else if(outcome == 0) // Adhered
          upper_bound = magnitude;
        magnitude = lower_bound + (upper_bound - lower_bound) / 2;
      }
      ++output_i;
    }
  }

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  std::cout << "time: " << timer << " ms " << (timer/1000) << " s" << std::endl;

	obj["grid"] = grid;
  std::ofstream File("pp.json");
  json::print(File,obj);
  std::ofstream File2("/home/marmaduke/pp.json");
  json::print(File2,obj);
}

void equillibriate(json::Object& obj)
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

  std::vector< value_type > times;
	vector_type prev(2*size+2);
	vector_type curr(2*size+2);

  times.push_back(9);
	times.push_back(10);

	bool swtch = false;
	int obs_count = 0;

  force_functor F(p);
  pp_observer O(prev,curr,obj,swtch,obs_count);

  value_type dt = .05;
	
	while(true)
	{
		integrate_times(make_controlled(p.abstol, p.reltol, stepper_type()),
                	F, init, times.begin(), times.end(), dt, O);
		
		vector_type diff(2*size+2);
		thrust::minus<value_type> op;
		thrust::transform(curr.begin(),curr.end(),prev.begin(),diff.begin(),op);
		value_type equil = sqrt(thrust::inner_product(
								diff.begin(),diff.begin()+2*size,diff.begin(),0.));
    value_type adhesion = sqrt(thrust::inner_product(
								diff.begin()+2*size,diff.end(),diff.begin()+2*size,0.));
    value_type d = sqrt(p.lambda*p.lambda + p.mu*p.mu);
		std::cout << "equil: " << equil << " adhesion: " << adhesion << " d: " << d << std::endl;
		if(abs(d - adhesion) < p.abstol or 
			(equil < p.abstol and adhesion < p.abstol))
		{	
			break;
		}
	}
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timer,start,stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  json::Array* temp = new json::Array;
  std::stringstream ss;
  ss << "tq" << obs_count;
  std::string temp2 = ss.str();
  for(int j = 0; j < 2*size+2; ++j)
	{
    json::Number* num = new json::Number;
    num->val = init[j];
    temp->push_back(num);
  }
  obj[temp2] = temp;

  std::cout << "time: " << timer << " ms " << (timer/1000) << " s" << std::endl;

  std::ofstream File("output.json");
  json::print(File,obj);
  std::ofstream File2("/home/marmaduke/output.json");
  json::print(File2,obj);
}

int main()
{
  json::Object obj = json::parse(std::cin,10);
  json::Number* num = as<json::Number>(obj["type"]);
	json::Number* device = as<json::Number>(obj["device"]);
	int t = num->val;
	int d = device->val;
	switch(t)
	{
		case 1:
			std::cout << "Pulloff Profile" << std::endl;
			pulloff_profile(obj);
			break;
		default:
			std::cout << "Equillibriate" << std::endl;
			equillibriate(obj);
	}
	return 0;
}
