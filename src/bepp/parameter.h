#include <cstdlib>

#include "defs.h"
#include "../json/json.h"

#ifndef PARAMETER_H
#define PARAMETER_H

struct parameter
{
  int m;
  int n;
  value_type beta;
  value_type len;
  value_type gamma;
  value_type epsilon;
  value_type epsilon_bottom;
	value_type epsilon_top;
	value_type sigma;
  value_type lambda;
  value_type mu;
	value_type pressure;

	bool have_delta, have_init;
	value_type* delta;
	value_type* init;
  
	value_type sub_h;
  value_type sub_x;
  value_type sub_y;
  int sub_count;

  value_type osub_h;
  value_type osub;
  int osub_count;

  value_type reltol;
  value_type abstol;

  parameter(json::Object& obj)
  {
    bool failed = false;

    if(obj.count("m") > 0)
    {
      json::Number* temp = as<json::Number>(obj["m"]);
      m = temp->val;
    }else
      failed = true;

    if(obj.count("n") > 0)
    {
      json::Number* temp = as<json::Number>(obj["n"]);
      n = temp->val;
    }else
      failed = true;

    if(obj.count("beta") > 0)
    {
      json::Number* temp = as<json::Number>(obj["beta"]);
      beta = temp->val;
    }else
      failed = true;

    if(obj.count("len") > 0)
    {
      json::Number* temp = as<json::Number>(obj["len"]);
      len = temp->val;
    }else
      failed = true;

    if(obj.count("gamma") > 0)
    {
      json::Number* temp = as<json::Number>(obj["gamma"]);
      gamma = temp->val;
    }else
      failed = true;

    if(obj.count("epsilon") > 0)
    {
      json::Number* temp = as<json::Number>(obj["epsilon"]);
      epsilon = temp->val;
    }else
      failed = true;

    if(obj.count("epsilon_bottom") > 0)
    {
      json::Number* temp = as<json::Number>(obj["epsilon_bottom"]);
      epsilon_bottom = temp->val;
    }else
      failed = true;
    
		if(obj.count("epsilon_top") > 0)
    {
      json::Number* temp = as<json::Number>(obj["epsilon_top"]);
      epsilon_top = temp->val;
    }else
      failed = true;
    
		if(obj.count("sigma") > 0)
    {
      json::Number* temp = as<json::Number>(obj["sigma"]);
      sigma = temp->val;
    }else
      failed = true;

    if(obj.count("lambda") > 0)
    {
      json::Number* temp = as<json::Number>(obj["lambda"]);
      lambda = temp->val;
    }else
      failed = true;

    if(obj.count("mu") > 0)
    {
      json::Number* temp = as<json::Number>(obj["mu"]);
      mu = temp->val;
    }else
      failed = true;

    if(obj.count("sub_h") > 0)
    {
      json::Number* temp = as<json::Number>(obj["sub_h"]);
      sub_h = temp->val;
    }else
      failed = true;

    if(obj.count("sub_x") > 0)
    {
      json::Number* temp = as<json::Number>(obj["sub_x"]);
      sub_x = temp->val;
    }else
      failed = true;

    if(obj.count("sub_y") > 0)
    {
      json::Number* temp = as<json::Number>(obj["sub_y"]);
      sub_y = temp->val;
    }else
      failed = true;

    if(obj.count("sub_count") > 0)
    {
      json::Number* temp = as<json::Number>(obj["sub_count"]);
      sub_count = temp->val;
    }else
      failed = true;

    if(obj.count("osub_h") > 0)
    {
      json::Number* temp = as<json::Number>(obj["osub_h"]);
      osub_h = temp->val;
    }else
      failed = true;

    if(obj.count("osub") > 0)
    {
      json::Number* temp = as<json::Number>(obj["osub"]);
      osub = temp->val;
    }else
      failed = true;

    if(obj.count("osub_count") > 0)
    {
      json::Number* temp = as<json::Number>(obj["osub_count"]);
      osub_count = temp->val;
    }else
      failed = true;

    if(obj.count("abstol") > 0)
    {
      json::Number* temp = as<json::Number>(obj["abstol"]);
      abstol = temp->val;
    }else
      failed = true;

    if(obj.count("reltol") > 0)
    {
      json::Number* temp = as<json::Number>(obj["reltol"]);
      reltol = temp->val;
    }else
      failed = true;
		
		if(obj.count("pressure") > 0)
		{
			json::Number* temp = as<json::Number>(obj["pressure"]);
			pressure = temp->val;
		}else
			failed = true;

		if(obj.count("delta") > 0)
    {
      json::Array* temp = as<json::Array>(obj["delta"]);
      delta = new value_type[temp->size()];
      for(int i = 0; i < temp->size(); ++i)
      {
        json::Number* temp2 = as<json::Number>((*temp)[i]);
        delta[i] = temp2->val;
      }
			have_delta = true;
    }else
			have_delta = false;

		if(obj.count("init") > 0)
		{
			json::Array* temp = as<json::Array>(obj["init"]);
			init = new value_type[temp->size()];
			for(int i = 0; i < temp->size(); ++i)
			{
				json::Number* temp2 = as<json::Number>((*temp)[i]);
				init[i] = temp2->val;
			}
			have_init = true;
		}else
			have_init = false;
    
		if(failed)
    {
      std::cout << "JSON file missing parameters" << std::endl;
      exit(1);
    }
  };

};

#endif

