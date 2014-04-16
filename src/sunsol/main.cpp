
#include <sstream>
#include <iostream>
#include <fstream>

#ifdef CUDA_TARGET
  //TODO
#else

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>             /* main integrator header file */
#include <cvode/cvode_dense.h>       /* use CVDENSE linear solver */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, */
#include <sundials/sundials_types.h> /* definition of realtype */
#include <sundials/sundials_math.h>  /* contains the macros ABS, SQR, and EXP */

#include "json/json.h"
#include "serial/force.h"
#include "defs.h"
#include "parameter.h"
#endif

#define ZERO RCONST(0)
#define ONE RCONST(1)
#define EIGHTH RCONST(.125)

static int equillibriate(parameter& params, json::Object& obj);
static N_Vector generate_init(parameter& params);
static void save(N_Vector y, json::Object& obj, uint count, uint nv_size);

int main()
{
  json::Object obj = json::parse(std::cin,10);
  json::Number num = *((json::Number*) obj["type"]);
  json::Number device = *((json::Number*) obj["device"]);
  int t = num.val;
  int d = device.val;

#ifdef CUDA_TARGET
  cudaSetDevice(d);
#endif

  parameter p(obj);
  switch(t)
  {
    default:
      equillibriate(p, obj);
  }
}

static void save(N_Vector y, json::Object& obj, uint count, uint nv_size)
{
  json::Array* a = new json::Array();
  for(uint i = 0; i < nv_size; ++i)
  {
    realtype x = NV_Ith_S(y, i);
    a->push_back(new json::Number(x));
  }
  std::stringstream ss;
  ss << "tq" << count;
  obj[ss.str()] = a;
}

static int equillibriate(parameter& params, json::Object& obj)
{
  realtype reltol = params.reltol, abstol = params.abstol;
  realtype movtol = params.movtol;
  realtype lambda = params.lambda, mu = params.mu;
  realtype t = ZERO, tout = EIGHTH, tpause = RCONST(100), equil;
  int flag, query;
  uint count = 0;
  N_Vector y = NULL, y_ref = NULL, v_equil = NULL;
  void *cvode_mem = NULL;

  uint particle_count = params.n * params.m + 1;
  uint nv_size = 2*particle_count;
  uint sub_index = 2 * params.n * params.m;

  force_wrapper Wrapper(params);

  if(params.have_init)
    y = N_VMake_Serial(nv_size, params.init);
  else
    y = generate_init(params);
  y_ref = N_VClone(y);
  v_equil = N_VClone(y);

  // CV_ADAMS, CV_FUNCTIONAL (Adams-Bashforth-Moulton)
  // CV_BDF, CV_NEWTON (Implicit Backward Difference Formula)
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  flag = CVodeInit(cvode_mem, force, ZERO, y);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
#ifdef CUDA_TARGET
  // TODO
#else
  /*  As of now I see no reason to use anything but the provided Dense solver in
      a serial setting. van der Waals removes any hope of a banded structure and
      any spareness will be too dynamic to be worth capturing.
      The Krylov iteration is however worth exploring.
  */
  flag = CVDense(cvode_mem, nv_size);
#endif

  flag = CVodeSetUserData(cvode_mem, &Wrapper);

  while(true)
  {
  #ifdef CUDA_TARGET
    // TODO, Might want to just generalize save function
  #else
    save(y, obj, count++, nv_size);
  #endif

    if(t >= tpause)
    {
      printf("tpause has been reached, end (0) or continue (1)?\n");
      scanf("%d",&query);
      if(query == 0)
        break;
      tpause += RCONST(100);
    }

    N_VAddConst(y, ZERO, y_ref); // Copy y to y_ref
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    tout += EIGHTH;

    N_VLinearSum(ONE, y, RCONST(-1), y_ref, v_equil); // Compute y - y_ref
    equil = N_VMaxNorm(v_equil);

  #ifdef CUDA_TARGET
    // TODO
  #else
    realtype sub_x = NV_Ith_S(v_equil, sub_index);
    realtype sub_y = NV_Ith_S(v_equil, sub_index+1);
  #endif

    // Scale by 8 because of jump in time (by an eighth)
    realtype adhesion = RCONST(8)*sqrt(sub_x*sub_x + sub_y*sub_y);
    realtype term_velocity = sqrt(lambda*lambda + mu*mu);

    if(abs(term_velocity- adhesion) < movtol or
      (equil < movtol and adhesion < movtol))
    {
      break;
    }
  }

  std::ofstream File("output.json");
  json::print(File,obj);

  N_VDestroy(y);
  N_VDestroy(y_ref);
  N_VDestroy(v_equil);
  CVodeFree(&cvode_mem);
}

static N_Vector generate_init(parameter& params)
{
  int n = params.n, m = params.m;
  uint size = n*m;
  uint nv_size = RCONST(2)*(n*m + ONE);
  N_Vector out = N_VNew_Serial(nv_size);

  realtype* temp_delta;

  if(not params.have_delta)
    temp_delta = new realtype[m];

  for(int j = 0; j < m; ++j)
  {
    if(params.have_delta)
      temp_delta[j] = params.delta[j];
    else
      temp_delta[j] = (realtype)j;
    for(int i = 0; i < n; ++i)
    {
      int index = i + n*j;
      NV_Ith_S(out, index) = temp_delta[j];
      NV_Ith_S(out, index + size) = (realtype)i + ONE;
    }
  }

  if(params.have_delta)
    delete temp_delta;
  else
  {
    params.delta = temp_delta;
  }

  return out;
}
