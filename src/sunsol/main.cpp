
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <cassert>
#include <ctime>

#ifdef CUDA_TARGET
  //TODO
#else

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <csignal>

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

json::Object* OBJ_PTR;

static int equillibriate(parameter& params, json::Object& obj);
static N_Vector generate_init(parameter& params);
static void save(N_Vector y, json::Object& obj, uint count, uint nv_size);
void handle_kill(int sig);

int main()
{
  json::Object obj = json::parse(std::cin,10);
  OBJ_PTR = &obj;
  json::Number num = *((json::Number*) obj["type"]);
  json::Number device = *((json::Number*) obj["device"]);
  int t = num.val;
  int d = device.val;

  std::signal(SIGINT, handle_kill);

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
  realtype t = ZERO, tout = RCONST(8), tsave = RCONST(1000), equil;
  int flag, query = 1;
  uint count = 0;
  N_Vector y = NULL, nhbd_ref = NULL, y_ref = NULL, v_equil = NULL;
  void *cvode_mem = NULL;

  uint particle_count = params.n * params.m + 1;
  uint nv_size = 2*particle_count;
  uint sub_index = 2 * params.n * params.m;

  if(params.have_init)
    y = N_VMake_Serial(nv_size, params.init);
  else
    y = generate_init(params);
  y_ref = N_VClone(y);
  v_equil = N_VClone(y);
  nhbd_ref = N_VClone(y);
  N_VAddConst(y, ZERO, nhbd_ref);

  force_wrapper Wrapper(params, nhbd_ref, y);

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
  flag = CVodeSetMaxNumSteps(cvode_mem, -1);
  //flag = CVodeSetMaxConvFails(cvode_mem, 1000);

  //force(RCONST(0), y, y_ref, &Wrapper);
  //N_VPrint_Serial(y_ref);
  //assert(false);

  clock_t start = clock();
  int counter = 0;
  while(counter < 10)
  {
    /*
    if(t >= tsave)
    {
      save(y, obj, count++, nv_size);
      tsave += RCONST(1000);
    } */

    N_VAddConst(y, ZERO, y_ref); // Copy y to y_ref
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    tout += RCONST(8);
    //std::cout << tout << " " << equil << " " << flag << std::endl;

    N_VLinearSum(ONE, y, RCONST(-1), y_ref, v_equil); // Compute y - y_ref

  #ifdef CUDA_TARGET
    // TODO
  #else
    realtype sub_x = NV_Ith_S(v_equil, sub_index);
    realtype sub_y = NV_Ith_S(v_equil, sub_index+1);
    NV_Ith_S(v_equil, sub_index) = 0;
    NV_Ith_S(v_equil, sub_index + 1) = 0;
  #endif

    equil = N_VMaxNorm(v_equil);

    // Scale by 8 because of jump in time (by an eighth)
    realtype adhesion = RCONST(.125)*sqrt(sub_x*sub_x + sub_y*sub_y);
    adhesion = params.sub_count == 0 ? ZERO : adhesion;
    realtype term_velocity = sqrt(lambda*lambda + mu*mu);

    std::cout << "tout: " << (tout - RCONST(8)) << " equil: " << equil << std::endl;

    if(abs(term_velocity - adhesion) < movtol or
      (equil < movtol and adhesion < movtol))
    {
      ++counter;
    }else
      counter = 0;
  }
  clock_t end = clock();
  std::cout << "Time: " << (end - start)/CLOCKS_PER_SEC << " seconds. "
    << "Simulation Time: " << (tout - RCONST(8)) << std::endl;

  N_VPrint_Serial(y);

  /*
  for(int i = 0, b = 0; i < Wrapper.nhbd_fiber.size(); ++i, b+=2)
  {
    if(Wrapper.mask[b] and not Wrapper.mask[b+1])
    {
      std::cout << "Upper:" << std::endl;
      std::cout << Wrapper.nhbd_fiber[i] << " " << Wrapper.nhbd_partner[i] << std::endl;
    }
    if(not Wrapper.mask[b] and Wrapper.mask[b+1])
    {
      std::cout << "Lower:" << std::endl;
      std::cout << Wrapper.nhbd_fiber[i] << " " << Wrapper.nhbd_partner[i] << std::endl;
    }
  }
  */

  save(y, obj, count++, nv_size);

  std::ofstream File("output.json");
  json::print(File,obj);

  N_VDestroy(y);
  N_VDestroy(y_ref);
  N_VDestroy(v_equil);
  N_VDestroy(nhbd_ref);
  CVodeFree(&cvode_mem);
}

static N_Vector generate_init(parameter& params)
{
  int n = params.n, m = params.m;
  uint size = n*m;
  uint nv_size = RCONST(2)*(n*m + ONE);
  N_Vector out = N_VNew_Serial(nv_size);

  realtype* temp_delta;

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

  NV_Ith_S(out, 2*size) = params.sub_x;
  NV_Ith_S(out, 2*size + 1) = params.sub_y;

  if(params.have_delta)
    delete temp_delta;
  else
  {
    params.delta = temp_delta;
  }

  return out;
}

void handle_kill(int sig)
{
  printf("Caught signal %d\nWriting partial JSON object data.\n",sig);
  std::ofstream File("terminated_output.json");
  json::print(File,*OBJ_PTR);
  exit(1);
}
