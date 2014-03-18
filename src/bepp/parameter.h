#include <cstdlib>

#include "defs.h"
#include "../json/json.h"

#ifndef PARAMETER_H
#define PARAMETER_H

struct parameter
{
  int m;
  int n;
  double beta;
  double len;
  double gamma;
  double epsilon;
  double sigma;
	double lambda;
	double mu;

  double sub_h;
  double sub_x;
	double sub_y;
	int sub_count;
  double * delta;

	double osub_h;
	double osub;
	int osub_count;

	double reltol;
	double abstol;

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
/*
    if(obj.count("delta") > 0)
    {
      json::Array* temp = as<json::Array>(obj["delta"]);
      cudaMalloc(&delta,sizeof(double)*temp->size());
      double* temp2 = new double[temp->size()];
      for(int i = 0; i < temp->size(); ++i)
      {
        json::Number* temp3 = as<json::Number>((*temp)[i]);
        temp2[i] = temp3->val;
      }
      cudaMemcpy(delta,temp2,sizeof(double)*temp->size(),cudaMemcpyHostToDevice);
    }
    */
		if(failed)
    {
      std::cout << "JSON file missing parameters" << std::endl;
      exit(1);
    }
  };

};

#endif

