#include <cstdlib>
#include <limits>

#include "defs.h"
#include "json/json.h"

#ifndef PARAMETER_H
#define PARAMETER_H

#define PROCESS_NUMBER_REQ(X, Y) \
      if(obj.count((X)) > 0) \
      { \
        json::Number* temp = as<json::Number>(obj[(X)]); \
        (Y) = temp->val; \
      }else \
        failed = true

#define PROCESS_NUMBER_OPT(X, Y, Z) \
      if(obj.count((X)) > 0) \
      { \
        json::Number* temp = as<json::Number>(obj[(X)]); \
        (Y) = temp->val; \
      }else \
        (Y) = (Z)

#define PROCESS_BOOL_OPT(X, Y, Z) \
      if(obj.count((X)) > 0) \
      { \
        json::Bool* temp = as<json::Bool>(obj[(X)]); \
        (Y) = temp->val; \
      }else \
        (Y) = (Z)

struct parameter
{
  int m;
  int n;
  realtype beta;
  realtype len;
  realtype gamma;
  realtype epsilon;
  realtype epsilon_bottom;
  realtype epsilon_top;
  realtype epsilon_subs;
  realtype sigma;
  realtype lambda;
  realtype mu;
  realtype pressure;

  bool f2f_switch;
  bool f2l_switch;
  bool f2u_switch;
  bool s2s_switch;

  realtype rcut;
  realtype rmax;

  bool have_delta, have_init;
  realtype* delta;
  realtype* init;

  realtype sub_h;
  realtype sub_x;
  realtype sub_y;
  int sub_count;

  realtype osub_h;
  realtype osub;
  int osub_count;

  realtype reltol;
  realtype abstol;
  realtype movtol;

  bool save_all;
  int save;
  realtype save_step;
  realtype tmax;

  parameter(json::Object& obj)
  {
    bool failed = false;

    PROCESS_NUMBER_REQ("m",m);
    PROCESS_NUMBER_REQ("n",n);
    PROCESS_NUMBER_REQ("sub_h",sub_h);
    PROCESS_NUMBER_REQ("sub_x",sub_x);
    PROCESS_NUMBER_REQ("sub_y",sub_y);
    PROCESS_NUMBER_REQ("sub_count",sub_count);
    PROCESS_NUMBER_REQ("osub",osub);
    PROCESS_NUMBER_REQ("osub_h",osub_h);
    PROCESS_NUMBER_REQ("osub_count",osub_count);
    PROCESS_NUMBER_REQ("rcut",rcut);
    PROCESS_NUMBER_REQ("rmax",rmax);

    PROCESS_NUMBER_OPT("beta",beta,ONE);
    PROCESS_NUMBER_OPT("len",len,ONE);
    PROCESS_NUMBER_OPT("gamma",gamma,ONE);
    PROCESS_NUMBER_OPT("epsilon",epsilon,ONE);
    PROCESS_NUMBER_OPT("epsilon_top",epsilon_top,ONE);
    PROCESS_NUMBER_OPT("epsilon_bottom",epsilon_bottom,ONE);
    PROCESS_NUMBER_OPT("epsilon_subs",epsilon_subs,ZERO);
    PROCESS_NUMBER_OPT("sigma",sigma,ONE);
    PROCESS_NUMBER_OPT("lambda",lambda,ZERO);
    PROCESS_NUMBER_OPT("mu",mu,ZERO);
    PROCESS_NUMBER_OPT("abstol",abstol,RCONST(1e-9));
    PROCESS_NUMBER_OPT("reltol",reltol,RCONST(1e-9));
    PROCESS_NUMBER_OPT("movtol",movtol,RCONST(1e-6));
    PROCESS_NUMBER_OPT("pressure",pressure,ZERO);
    PROCESS_NUMBER_OPT("save_step",save_step,RCONST(.1));
    PROCESS_NUMBER_OPT("tmax",tmax,std::numeric_limits<realtype>::max());
    PROCESS_NUMBER_OPT("save",save,0);

    PROCESS_BOOL_OPT("f2f_switch",f2f_switch,true);
    PROCESS_BOOL_OPT("f2u_switch",f2u_switch,true);
    PROCESS_BOOL_OPT("f2l_switch",f2l_switch,true);
    PROCESS_BOOL_OPT("s2s_switch",s2s_switch,false);
    PROCESS_BOOL_OPT("save_all",save_all,true);

    if(obj.count("delta") > 0)
    {
      json::Array* temp = as<json::Array>(obj["delta"]);
      delta = new realtype[temp->size()];
      for(uint i = 0; i < temp->size(); ++i)
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
      init = new realtype[temp->size()];
      for(uint i = 0; i < temp->size(); ++i)
      {
        json::Number* temp2 = as<json::Number>((*temp)[i]);
        init[i] = temp2->val;
      }
      have_init = true;
    }else
      have_init = false;

    if(have_init and not have_delta)
      failed = true;

    if(failed)
    {
      std::cout << "JSON file missing parameters" << std::endl;
      exit(1);
    }
  };

};

#endif

