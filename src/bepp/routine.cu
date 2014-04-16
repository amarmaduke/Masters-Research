#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <cmath>

#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include <boost/ref.hpp>

#include "force.h"
#include "defs.h"
#include "../json/json.h"

using namespace boost::numeric::odeint;

template<typename T>
struct absolute_value : public thrust::unary_function<T,T>
{
  __host__ __device__
  T operator()(const T& x) const
  {
    return x < T(0) ? -x : x;
  }
};

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

struct pp_observer_p
{
  vector_type& prev;

  pp_observer_p(vector_type& p) : prev(p) { }

  template<typename State >
  void operator() (const State&x, value_type t)
  {
    if(t == 9)
    {
      thrust::copy(x.begin(),x.end(),prev.begin());
    }
  }
};

void pulloff_profile_p(parameter& p, json::Object& obj, vector_type& init)
{
  int size = p.m*p.n;
  int total_size = 2*size+2;
  cudaEvent_t start, stop;
  float timer;

  json::Number* m_tol = as<json::Number>(obj["movtol"]);
  double movtol = m_tol->val;

  std::vector< std::vector<value_type> > theta_ranges;
  json::Array* range = as<json::Array>(obj["range"]);
  for(int i = 0; i < SIM_COUNT; ++i)
  {
    std::vector<value_type> temp;
    json::Array* r_temp = as<json::Array>((*range)[i]);
    for(int j = 0; j < 3; ++j)
    {
      json::Number* n_temp = as<json::Number>((*r_temp)[j]);
      temp.push_back(n_temp->val);
    }
    theta_ranges.push_back(temp);
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);

  std::vector<value_type> times;
  vector_type prev(SIM_COUNT*total_size);
  times.push_back(9); times.push_back(10);

  force_functor2 F(p);
  pp_observer_p O(prev);
  vector_type V = vector_type(SIM_COUNT*total_size);
  for(int i = 0; i < SIM_COUNT; ++i)
  {
    thrust::copy(init.begin(),init.end(),V.begin()+i*total_size);
  }

  json::Array* grid = new json::Array();

  double linear_step[SIM_COUNT];
  double magnitude[SIM_COUNT];
  double lower_bound[SIM_COUNT];
  double upper_bound[SIM_COUNT];
  bool have_upper[SIM_COUNT];
  bool have_lower[SIM_COUNT];
  double theta_begin[SIM_COUNT];
  double theta_end[SIM_COUNT];
  double theta_h[SIM_COUNT];
  double t[SIM_COUNT];
  int outcome[SIM_COUNT];

  for(int i = 0; i < SIM_COUNT; ++i)
  {
    linear_step[i] = 250;
    magnitude[i] = 250;
    lower_bound[i] = 0;
    upper_bound[i] = 100;
    have_lower[i] = true;
    have_upper[i] = false;
    theta_begin[i] = theta_ranges[i][0];
    theta_end[i] = theta_ranges[i][2];
    theta_h[i] = theta_ranges[i][1];
    t[i] = theta_begin[i];
    outcome[i] = -1;
  }

  double f_t = 0, time = 0, p_time = 0;
  double dt = 1e-4, p_dt = 0;
  bool SIM_COMPLETE = false;

  double l[SIM_COUNT], m[SIM_COUNT];
	double adh_t[SIM_COUNT] = {0};
  bool found[SIM_COUNT] = {false};
  bool done[SIM_COUNT] = {false};
  bool complete[SIM_COUNT] = {false};

  controlled_runge_kutta< stepper_type > stepper
          = make_controlled(p.abstol,p.reltol,stepper_type());

  while(not SIM_COMPLETE)
  {
    SIM_COMPLETE = true;
    for(int i = 0; i < SIM_COUNT; ++i)
    {
      if(t[i] > theta_end[i])
        complete[i] = true;
      SIM_COMPLETE = SIM_COMPLETE and complete[i];
      if(found[i])
      {
        t[i] += theta_h[i];
        found[i] = false;
        have_lower[i] = true;
        lower_bound[i] = 0;
        have_upper[i] = false;
      }

      l[i] = -magnitude[i]*sin(PI*t[i]/180);
      m[i] = magnitude[i]*cos(PI*t[i]/180);
      if(complete[i])
			{
				l[i] = 0;
				m[i] = 0;
			}
			if(t[i] != theta_begin[i])
        linear_step[i] = 50;
      F.lambda[i] = l[i];
      F.mu[i] = m[i];
    }

    if(SIM_COMPLETE)
      break;


    //integrate_times(make_controlled(p.abstol, p.reltol, stepper_type()),
    //    F, V, times.begin(), times.end(), dt, O);
    controlled_step_result result = success;
    thrust::copy(V.begin(),V.end(),prev.begin());
    p_time = time;
    while(time - p_time < .125)
    {
      do
      {
        if(time - p_time + 1.1*dt > .125)
        {
          dt = .125 - time + p_time;
        }
        p_dt = dt;
        result = stepper.try_step(F,V,f_t,dt);
        if(dt < 1e-30)
        {
          const char * error_string = "dt is too small.";
          throw std::overflow_error( error_string );
        }
      }while(result != success);
      time += p_dt;
    }
    time = 0;

    for(int i = 0; i < SIM_COUNT; ++i)
    {
      vector_type diff(total_size);
      thrust::minus<value_type> op;
      thrust::transform(V.begin()+i*total_size,
                        V.begin()+(i+1)*total_size,
                        prev.begin()+i*total_size,
                        diff.begin(), op);

      value_type equil = thrust::transform_reduce(
                diff.begin(),diff.begin()+2*size,absolute_value<value_type>(),
                0.0,thrust::maximum<value_type>());

      value_type adhesion = sqrt(thrust::inner_product(
                  diff.begin()+2*size,diff.end(),diff.begin()+2*size,0.0));
			
			value_type adh_r = sqrt(thrust::inner_product(
									V.begin()+2*size+i*total_size,
									V.begin()+2*size+2+i*total_size,
									V.begin()+2*size+i*total_size,0.0));

      adhesion = 8.0*adhesion;

      value_type d = sqrt(l[i]*l[i] + m[i]*m[i]);

			if(adhesion < movtol or adhesion/adh_r < movtol)
				adh_t[i] += .125;
			else
				adh_t[i] = 0;

			if(not complete[i])
				std::cout << i << ": theta: " << t[i] << " " << equil << " " << (adhesion-d) << std::endl;
      else
				std::cout << i << ": complete." << std::endl;
			if(abs(adhesion - d) < movtol and not complete[i])
      { // Pulled off
        done[i] = true;
        json::Array* temp = new json::Array();
        temp->push_back(new json::Number(t[i]));
        temp->push_back(new json::Number(l[i]));
        temp->push_back(new json::Number(m[i]));
        temp->push_back(new json::Number(0));
        grid->push_back(temp);
        outcome[i] = 0;
        std::cout << "Pulled off. t: " << t[i] << " m: " << m[i] << " l: " << l[i] << std::endl;
      }else if(((equil < movtol and adhesion < movtol) or (adh_t[i] >= 1)) and not complete[i])
      { // Adhered
        done[i] = true;
        json::Array* temp = new json::Array();
        temp->push_back(new json::Number(t[i]));
        temp->push_back(new json::Number(l[i]));
        temp->push_back(new json::Number(m[i]));
        temp->push_back(new json::Number(1));
        grid->push_back(temp);
        outcome[i] = 1;
        std::cout << "Adhered. t: " << t[i] << " m: " << m[i] << " l: " << l[i] << std::endl;
      }

      if(done[i])
      {
        thrust::copy(init.begin(),init.end(),V.begin()+i*total_size);
        if(have_upper[i] and have_lower[i] and abs(upper_bound[i] - lower_bound[i]) < .5)
        {
          found[i] = true;
        }

        if(not have_lower[i] or not have_upper[i])
        {
          if(outcome[i] == 1) // Adhered
          {
            lower_bound[i] = magnitude[i];
            have_lower[i] = true;
            if(not have_upper[i])
              magnitude[i] = max(magnitude[i] + linear_step[i],0.);
          }else if(outcome[i] == 0) // Pulled Off
          {
            upper_bound[i] = magnitude[i];
            have_upper[i] = true;
            if(not have_lower[i])
              magnitude[i] = max(magnitude[i] - linear_step[i],0.);
          }
        }

        if(have_lower[i] and have_upper[i])
        {
          if(outcome[i] == 1) // Adhered
            lower_bound[i] = magnitude[i];
          else if(outcome[i] == 0) // Pulled Off
            upper_bound[i] = magnitude[i];
          magnitude[i] = lower_bound[i] + (upper_bound[i] - lower_bound[i]) / 2;
        }
        done[i] = false;
				adh_t[i] = 0;
      }
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

void pulloff_profile(parameter& p, json::Object& obj, vector_type& init)
{
  int size = p.m*p.n;
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

  json::Number* m_tol = as<json::Number>(obj["movtol"]);
  double movtol = m_tol->val;

  double linear_step = 50;
  int output_i = 1;
  double magnitude = 450;
  double lower_bound = 0;
  double upper_bound = 100;
  bool have_upper = false, have_lower = false;

  double dt = .001, f_t = 0, time = 0, p_time = 0, adh_t = 0;
  double theta_begin = 40, theta_end = 45, theta_h = 5;

  controlled_runge_kutta< stepper_type > stepper
          = make_controlled(p.abstol,p.reltol,stepper_type());

  for(int t = theta_begin; t <= theta_end; t += theta_h)
  {
    bool found = false;
    have_lower = have_upper = false;

    if(t != theta_begin)
    {
      linear_step = 5;
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

			adh_t = 0;
    	controlled_step_result result = success;
      while(not done)
      {
				
    		double p_dt;
    		while(time - p_time < .125)
   			{
      		do
     	 		{
        		// Time is irrelevant for us.
        		if(time - p_time + 1.1*dt > .125)
       	 		{
          		dt = .125 - time + p_time;
        		}
        		p_dt = dt;
        		result = stepper.try_step(F,V,f_t,dt);
        		if(dt < 1e-25)
        		{
          		const char * error_string = "dt is too small.";
          		throw std::overflow_error( error_string );
        		}
      		}while(result != success);
      		time += p_dt;
    		}
        
      	vector_type diff(2*size+2);
      	thrust::minus<value_type> op;
      	thrust::transform(V.begin(),V.end(),prev.begin(),diff.begin(),op);

      	value_type equil = thrust::transform_reduce(
                  diff.begin(),diff.begin()+2*size,absolute_value<value_type>(),
                  0.0,thrust::maximum<value_type>());

      	value_type adhesion = sqrt(thrust::inner_product(
                  diff.begin()+2*size,diff.end(),diff.begin()+2*size,0.0));

				value_type adh_r = sqrt(thrust::inner_product(
									V.begin()+2*size,V.end(),V.begin()+2*size,0.0));

      	adhesion = adhesion*8;

				if(adhesion < movtol or adhesion/adh_r < movtol)
					adh_t += .125;
				else
					adh_t = 0;

      	value_type d = sqrt(l*l + m*m);

        if(abs(adhesion - d) < movtol)
        { // Pulled off
          done = true;
          json::Array* temp = new json::Array();
          temp->push_back(new json::Number(t));
          temp->push_back(new json::Number(l));
          temp->push_back(new json::Number(m));
          temp->push_back(new json::Number(0));
          grid->push_back(temp);
          outcome = 0;
          std::cout << "Pulled off. l: " << l << " m: " << m << std::endl;
        }else if((equil < movtol and adhesion < movtol) or adh_t >= 1)
        { // Adhered
          done = true;
          json::Array* temp = new json::Array();
          temp->push_back(new json::Number(t));
          temp->push_back(new json::Number(l));
          temp->push_back(new json::Number(m));
          temp->push_back(new json::Number(1));
          grid->push_back(temp);
          outcome = 1;
          std::cout << "Adhered. l: " << l << " m: " << m << std::endl;
        }
      }
      thrust::copy(init.begin(),init.end(),V.begin());

      if(have_upper and have_lower and abs(upper_bound - lower_bound) < .5)
      {
        found = true;
      }

      if(not have_lower or not have_upper)
      {
        if(outcome == 1) // Adhered
        {
          lower_bound = magnitude;
          have_lower = true;
          if(not have_upper)
            magnitude = magnitude + linear_step;
        }else if(outcome == 0) // Pulled Off
        {
          upper_bound = magnitude;
          have_upper = true;
          if(not have_lower)
            magnitude = magnitude - linear_step;
        }
      }

      if(have_lower and have_upper)
      {
        if(outcome == 1) // Adhered
          lower_bound = magnitude;
        else if(outcome == 0) // Pulled Off
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
  std::ofstream File("old_pp.json");
  json::print(File,obj);
  std::ofstream File2("/home/marmaduke/old_pp.json");
  json::print(File2,obj);
}

void equillibriate(parameter& p, json::Object& obj, vector_type& init)
{
  int size = p.m*p.n;
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

  json::Number* m_tol = as<json::Number>(obj["movtol"]);
  double movtol = m_tol->val;


  force_functor F(p);
  pp_observer O(prev,curr,obj,swtch,obs_count);

  controlled_runge_kutta< stepper_type > stepper
          = make_controlled(p.abstol,p.reltol,stepper_type());

  value_type dt = 1e-6;
  double time = 0, p_time = 0, s_time = 0, f_t = 0;
  double adh_t = 0;
	int count = 0;

  thrust::copy(init.begin(),init.end(),prev.begin());

  while(true)
  {
    controlled_step_result result = success;

    p_time = time;
    thrust::copy(init.begin(),init.end(),prev.begin());
    double p_dt;
    while(time - p_time < .125)
    {
      do
      {
        // Time is irrelevant for us.
        if(time - p_time + 1.1*dt > .125)
        {
          dt = .125 - time + p_time;
        }
        p_dt = dt;
        result = stepper.try_step(F,init,f_t,dt);
        if(dt < 1e-25)
        {
          const char * error_string = "dt is too small.";
          throw std::overflow_error( error_string );
        }
      }while(result != success);
      time += p_dt;
    }
    printf("%3.4f, %3.4f\n",(time - p_time), time);

    if(time - s_time > .1)
    {
      std::cout << "s_time: " << s_time << " time: " << time << std::endl;
      json::Array* a = new json::Array();
      for(int i = 0; i < init.size(); ++i)
      {
        double x = init[i];
        a->push_back(new json::Number(x));
      }
      std::stringstream ss;
      ss << "tq" << count;
      obj[ss.str()] = a;
      ++count;
      s_time = time;
    }
      vector_type diff(2*size+2);
      thrust::minus<value_type> op;
      thrust::transform(init.begin(),init.end(),prev.begin(),diff.begin(),op);

      value_type equil = thrust::transform_reduce(
                  diff.begin(),diff.begin()+2*size,absolute_value<value_type>(),
                  0.0,thrust::maximum<value_type>());

      value_type equil_r = thrust::transform_reduce(
                  init.begin(),init.begin()+2*size,absolute_value<value_type>(),
                  0.0,thrust::maximum<value_type>());
      
			value_type adhesion = sqrt(thrust::inner_product(
                  diff.begin()+2*size,diff.end(),diff.begin()+2*size,0.0));

			value_type adhesion_r = sqrt(thrust::inner_product(
                  init.begin()+2*size,init.end(),init.begin()+2*size,0.0));

      value_type d = sqrt(p.lambda*p.lambda + p.mu*p.mu);

			if(adhesion < movtol or adhesion/adhesion_r < movtol)
				adh_t += .125;
			else
				adh_t = 0;

      std::cout << "abs_eq: " << equil << " rel_eq: " << (equil/equil_r) << " abs_adh " << adhesion << " rel_adh: " << (adhesion/adhesion_r) << " d: " << d << " adh_t: " << adh_t << std::endl;

      if(abs(d - adhesion) < movtol or (adhesion < movtol and equil < movtol) or (adhesion/adhesion_r < movtol and equil/equil_r < movtol) or time > 10 or adh_t >= 1)
      {
        break;
      }
      p_time = time;
      thrust::copy(init.begin(),init.end(),prev.begin());
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

  cudaSetDevice(d);

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
        init[index] = d_delta[j] - ((double)i);
        init[index+size] = 1;
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

  switch(t)
  {
    case 1:
      std::cout << "Pulloff Profile" << std::endl;
      pulloff_profile(p,obj,init);
      break;
    case 2:
      std::cout << "Pulloff Profile Parallel" << std::endl;
      pulloff_profile_p(p,obj,init);
      break;
    default:
      std::cout << "Equillibriate" << std::endl;
      equillibriate(p,obj,init);
  }
  return 0;
}
