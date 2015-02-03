
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
json::Array* GRID_PTR;
bool save_all = true;

int pulloff_profile(parameter& params, json::Object& obj);
int pulloff_grid(parameter& params, json::Object& obj);
int pushon_grid(parameter& params, json::Object& obj);
int free_standing(parameter& params, json::Object& obj);
int pulloff_adh_bias(parameter& params, json::Object &obj);
int equillibriate(parameter& params, json::Object& obj, int sindex,
                  realtype scount, uint& bot_count, uint& top_count);
int equillibriate_fixed(parameter& params, json::Object& obj, int sindex,
                  realtype scount, uint& bot_count, uint& top_count);
N_Vector generate_init(parameter& params);
N_Vector copy_init(parameter& params);
realtype round_nearest_multiple(realtype t, realtype multiple);
void save(N_Vector y, json::Object& obj, std::string key, uint nv_size);
void save_grid( json::Array* grid, realtype t, realtype l, realtype m,
                int out, int sim, uint bot_count, uint top_count);
void handle_kill(int sig);
void check_CVode_error(int flag);
void compute_adhesion(N_Vector v, parameter& p,uint& bot_count,uint& top_count);

int main()
{
  json::Object obj = json::parse(std::cin,10);
  OBJ_PTR = &obj;
  json::Number num = *((json::Number*) obj["type"]);
  json::Number device = *((json::Number*) obj["device"]);
  json::Number save_num = *((json::Number*) obj["save"]);
  int t = num.val;
  //int d = device.val;

  std::signal(SIGINT, handle_kill);

#ifdef CUDA_TARGET
  cudaSetDevice(d);
#endif

  parameter p(obj);
  save_all = p.save_all;
  realtype save_step = p.save_step;
  clock_t start = clock();
  switch(t)
  {
    case 1:
      pulloff_profile(p, obj);
      break;
    case 2:
      pushon_grid(p, obj);
      break;
    case 3:
      pulloff_adh_bias(p, obj);
      break;
    case 4:
      free_standing(p, obj);
      break;
    default:
      uint b = 0, t = 0;
      equillibriate_fixed(p, obj, 0, save_step, b, t);
  }
  clock_t end = clock();
  std::cout << "Execution time: " << (end - start)/CLOCKS_PER_SEC
            << " seconds." << std::endl;

  std::stringstream ss;
  ss << save_num.val;
  std::string file_name(ss.str() + ".json");
  char* char_name = &file_name[0];
  std::ofstream File(char_name);
  json::print(File,obj);
}

void save(N_Vector y, json::Object& obj, std::string key, uint nv_size)
{
  if (not save_all)
    return;
  json::Array* a = new json::Array();
  for(uint i = 0; i < nv_size; ++i)
  {
    realtype x = NV_Ith_S(y, i);
    a->push_back(new json::Number(x));
  }
  obj[key] = a;
}

void save_grid( json::Array* grid, realtype t, realtype l, realtype m,
                int out, int sim, uint bot_count, uint top_count)
{
  json::Array* temp = new json::Array();
  temp->push_back(new json::Number(t));
  temp->push_back(new json::Number(l));
  temp->push_back(new json::Number(m));
  temp->push_back(new json::Number(out));
  temp->push_back(new json::Number(sim));
  temp->push_back(new json::Number(bot_count));
  temp->push_back(new json::Number(top_count));
  grid->push_back(temp);
}

int free_standing(parameter& params, json::Object& obj)
{
  json::Array* grid = new json::Array();
  GRID_PTR = grid;
  int sim_index = 0, outcome;
  uint bcount = 0, tcount = 0;

  std::vector< realtype > r1;
  json::Array* range = as<json::Array>(obj["range"]);
  for(uint i = 0; i < range->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range)[i]);
    realtype temp = numtemp->val;
    r1.push_back(temp);
  }

  std::vector< realtype > r2;
  json::Array* range2 = as<json::Array>(obj["range2"]);
  for(uint i = 0; i < range2->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range2)[i]);
    realtype temp = numtemp->val;
    r2.push_back(temp);
  }

  if(r2.size() == 0)
  {
    r2.push_back(1);
    r2.push_back(1);
    r2.push_back(10);
  }

  for(realtype i = r1[0]; i < r1[2]; i += r1[1])
  {
    std::cout << "Trying... i = " << i << std::endl;
    for(realtype j = r2[0]; j < r2[2]; j += r2[1])
    {
      std::cout << "\tTrying... j = " << j << std::endl;

      params.beta = i;
      params.epsilon_bottom = j;
      bcount = 0; tcount = 0;
      outcome = equillibriate(params, obj, sim_index, 0, bcount, tcount);

      save_grid(grid, i, j, 0, outcome, sim_index++, bcount, tcount);
    }
  }

  obj["grid"] = grid;
  return 0;
}

int pushon_grid(parameter& params, json::Object& obj)
{
  json::Array* grid = new json::Array();
  GRID_PTR = grid;
  int sim_index = 0, outcome;
  uint bcount = 0, tcount = 0;

  std::vector< realtype > tr;
  json::Array* range = as<json::Array>(obj["range"]);
  for(uint i = 0; i < range->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range)[i]);
    realtype temp = numtemp->val;
    tr.push_back(temp);
  }

  std::vector< realtype > mr;
  json::Array* range2 = as<json::Array>(obj["range2"]);
  for(uint i = 0; i < range2->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range2)[i]);
    realtype temp = numtemp->val;
    mr.push_back(temp);
  }

  if(mr.size() == 0)
  {
    mr.push_back(1);
    mr.push_back(1);
    mr.push_back(10);
  }

  for(realtype theta = tr[0]; theta < tr[2]; theta += tr[1])
  {
    std::cout << "Trying... theta = " << theta << std::endl;
    for(realtype magnitude = mr[0]; magnitude < mr[2]; magnitude += mr[1])
    {
      std::cout << "Trying... magnitude = " << magnitude << std::endl;
      realtype l = magnitude*sin(PI*theta/RCONST(180));
      realtype m = magnitude*cos(PI*theta/RCONST(180));

      params.lambda = l;
      params.mu = m;
      bcount = 0; tcount = 0;
      outcome = equillibriate(params, obj, sim_index, 0, bcount, tcount);

      save_grid(grid, theta, l, m, outcome, sim_index++, bcount, tcount);
    }
  }

  obj["grid"] = grid;
  return 0;
}

int pulloff_adh_bias(parameter& params, json::Object &obj)
{
  std::vector< realtype > theta_ranges;
  json::Array* range = as<json::Array>(obj["range"]);
  for(uint i = 0; i < range->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range)[i]);
    realtype temp = numtemp->val;
    theta_ranges.push_back(temp);
  }

  json::Array* grid = new json::Array();
  GRID_PTR = grid;

  bool have_upper = false, have_lower = false;
  realtype upper_bound, lower_bound = 0, magnitude = 10;
  uint tcount = 0, bcount = 0;
  int sim_index = 0;

  for(uint k = 0; k < theta_ranges.size(); k+=3)
  {
    realtype theta = theta_ranges[k]; // theta begin
    realtype theta_h = theta_ranges[k+1];
    realtype theta_end = theta_ranges[k+2];
    while(theta <= theta_end)
    {
      bool found = false;
      have_upper = false;
      have_lower = false;
      int outcome;
      std::cout << "Searching... theta = " << theta << std::endl;
      magnitude = 10;

      while(not found)
      {
        std::cout << "Trying... magnitude = " << magnitude << std::endl;
        realtype l = -magnitude*sin(PI*theta/RCONST(180));
        realtype m = magnitude*cos(PI*theta/RCONST(180));

        params.lambda = l;
        params.mu = m;

        outcome = equillibriate(params, obj, sim_index, 0, bcount, tcount);
        uint t_bcount = bcount, t_tcount = tcount;

        if(outcome == 0)
        {
          realtype mag1 = magnitude + RCONST(.25);
          int test1;

          std::cout << "Trying... magnitude = " << mag1 << std::endl;
          realtype l1 = -mag1*sin(PI*theta/RCONST(180));
          realtype m1 = mag1*cos(PI*theta/RCONST(180));
          params.lambda = l1;
          params.mu = m1;
          test1 = equillibriate(params, obj, sim_index, 0, bcount, tcount);

          if(test1 == 1)
          {
            std::cout << "Isolated Pulloff at " << magnitude << std::endl;
            save_grid(grid, theta, l1, m1, test1, sim_index++, bcount, tcount);
            bcount = 0; tcount = 0;
            magnitude = mag1; l = l1; m = m1; outcome = test1;
          }else if(test1 == 0)
          {
            save_grid(grid, theta, l, m, 2, sim_index++, bcount, tcount);
            bcount = 0; tcount = 0;
          }
        }

        save_grid(grid, theta, l, m, outcome, sim_index++, t_bcount, t_tcount);
        bcount = 0; tcount = 0;

        if(have_upper and have_lower and std::abs(upper_bound - lower_bound) < .1)
        {
          found = true;
        }

        if(not have_upper or not have_lower)
        {
          if(outcome == 1) // Adhered
          {
            lower_bound = magnitude;
            have_lower = true;
            if(not have_upper)
              magnitude += 10;
          }else if(outcome == 0) // Pulled off
          {
            upper_bound = magnitude;
            have_upper = true;
            if(not have_lower)
              magnitude /= 1.1;
          }
        }

        if(have_upper and have_lower)
        {
          if(outcome == 1) // Adhered
          {
            lower_bound = magnitude;
          }else if(outcome == 0) // Pulled off
          {
            upper_bound = magnitude;
          }
          magnitude = lower_bound + (upper_bound - lower_bound) / TWO;
        }
      }
      theta += theta_h;
    }
  }

  obj["grid"] = grid;
  return 0;
}

int pulloff_profile(parameter& params, json::Object& obj)
{

  std::vector< realtype > theta_ranges;
  json::Array* range = as<json::Array>(obj["range"]);
  for(uint i = 0; i < range->size(); ++i)
  {
    json::Number* numtemp = as<json::Number>((*range)[i]);
    realtype temp = numtemp->val;
    theta_ranges.push_back(temp);
  }

  json::Array* grid = new json::Array();
  GRID_PTR = grid;

  bool have_upper = false, have_lower = false, initial_guess = true;
  realtype upper_bound, lower_bound = 0, magnitude = 50;
  uint tcount = 0, bcount = 0;
  int sim_index = 0;

  for(uint k = 0; k < theta_ranges.size(); k+=3)
  {
    realtype theta = theta_ranges[k]; // theta begin
    realtype theta_h = theta_ranges[k+1];
    realtype theta_end = theta_ranges[k+2];
    while(theta <= theta_end)
    {
      bool found = false;
      have_upper = false;
      have_lower = false;
      int outcome;
      std::cout << "Searching... theta = " << theta << std::endl;

      while(not found)
      {
        std::cout << "Trying... magnitude = " << magnitude << std::endl;
        realtype l = -magnitude*sin(PI*theta/RCONST(180));
        realtype m = magnitude*cos(PI*theta/RCONST(180));

        params.lambda = l;
        params.mu = m;

        bcount = 0; tcount = 0;
        outcome = equillibriate(params, obj, sim_index, 0, bcount, tcount);
        uint t_bcount = bcount, t_tcount = tcount;

        if(outcome == 0)
        {
          realtype mag1 = magnitude + RCONST(.25);
          int test1;

          std::cout << "Trying... magnitude = " << mag1 << std::endl;
          realtype l1 = -mag1*sin(PI*theta/RCONST(180));
          realtype m1 = mag1*cos(PI*theta/RCONST(180));
          params.lambda = l1;
          params.mu = m1;

          bcount = 0; tcount = 0;
          test1 = equillibriate(params, obj, sim_index, 0, bcount, tcount);

          if(test1 == 1)
          {
            std::cout << "Isolated Pulloff at " << magnitude << std::endl;
            save_grid(grid, theta, l, m, 2, sim_index++, t_bcount, t_tcount);
            t_bcount = bcount; t_tcount = tcount;
            magnitude = mag1; l = l1; m = m1; outcome = test1;
          }else if(test1 == 0)
          {
            save_grid(grid, theta, l1, m1, test1, sim_index++, bcount, tcount);
          }
        }

        save_grid(grid, theta, l, m, outcome, sim_index++, t_bcount, t_tcount);

        if(have_upper and have_lower and std::abs(upper_bound - lower_bound)
                                      /std::min(upper_bound, lower_bound) < .01)
        {
          found = true;
        }

        if(not have_upper or not have_lower)
        {
          if(outcome == 1) // Adhered
          {
            lower_bound = magnitude;
            have_lower = true;
            if(not have_upper and not initial_guess)
              magnitude *= 1.1;
            if(not have_upper and initial_guess)
              magnitude *= 2;
          }else if(outcome == 0) // Pulled off
          {
            upper_bound = magnitude;
            have_upper = true;
            if(not have_lower and not initial_guess)
              magnitude /= 1.1;
            if(not have_lower and initial_guess)
              magnitude /= 2;
          }
        }

        if(have_upper and have_lower)
        {
          if(outcome == 1) // Adhered
          {
            lower_bound = magnitude;
          }else if(outcome == 0) // Pulled off
          {
            upper_bound = magnitude;
          }
          magnitude = lower_bound + (upper_bound - lower_bound) / TWO;
        }
      }
      theta += theta_h;
      initial_guess = false;
    }
  }

  obj["grid"] = grid;
  return 0;
}

int equillibriate(parameter& params, json::Object& obj, int sindex,
                  realtype scount, uint& bot_count, uint& top_count)
{
  realtype reltol = params.reltol, abstol = params.abstol;
  realtype movtol = params.movtol;
  realtype lambda = params.lambda, mu = params.mu;
  //realtype magnitude = sqrt(lambda*lambda + mu*mu);
  //realtype hmax = magnitude / (params.rmax - params.rcut) / RCONST(1.0001);
  realtype tp = ZERO, t = ZERO, tout = RCONST(8), equil;
  int flag;
  uint count = 0;
  N_Vector y = NULL, nhbd_ref = NULL, y_ref = NULL, nv_temp = NULL;
  void *cvode_mem = NULL;

  uint particle_count = params.n * params.m + 1;
  uint nv_size = 2*particle_count;
  uint sub_index = 2 * params.n * params.m;
  uint size = params.n * params.m;

  if(params.have_init)
    y = copy_init(params);
  else
    y = generate_init(params);
  y_ref = N_VClone(y);
  nv_temp = N_VClone(y);
  nhbd_ref = N_VClone(y);
  N_VAddConst(y, ZERO, nhbd_ref);

  force_wrapper Wrapper(params, y);

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
  flag = CVodeSetMaxConvFails(cvode_mem, 100);
//  flag = CVodeSetMaxStep(cvode_mem, hmax);

  //force(RCONST(0), y, y_ref, &Wrapper);
  //N_VPrint_Serial(y_ref);
  //assert(false);

  if(scount > 0)
  {
    std::stringstream ss;
    ss << "sindex" << sindex << "tq" << count++;
    save(y, obj, ss.str(), nv_size);
  }

  //clock_t start = clock();
  int counter = 0, outcome = 2;
  while(counter < 10)
  {

    N_VAddConst(y, ZERO, y_ref);
    tp = t;
    while(t < tout)
    {
      flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);
      if(flag < 0)
      {
        check_CVode_error(flag);
        outcome = -1;
        N_VPrint_Serial(y);
        exit(-1);
      }

      N_VLinearSum(1, y, -1, nhbd_ref, nv_temp);
      N_VProd(nv_temp, nv_temp, nv_temp); // Square componentwise
      realtype max_dist = ZERO, nmax_dist = ZERO;
      for(int i = 0; i < size; ++i)
      {
        realtype tx = NV_Ith_S(nv_temp, i);
        realtype ty = NV_Ith_S(nv_temp, i + size);
        realtype dist = sqrt(tx*tx + ty*ty);
        if(dist >= max_dist)
        {
          nmax_dist = max_dist;
          max_dist = dist;
        }
      }
      // Also check the upper substrate point.
      realtype tx = NV_Ith_S(nv_temp, 2*size);
      realtype ty = NV_Ith_S(nv_temp, 2*size + 1);
      realtype dist = sqrt(tx*tx + ty*ty);
      if(dist >= max_dist)
      {
        nmax_dist = max_dist;
        max_dist = dist;
      }

      if(max_dist + nmax_dist >= params.rmax - params.rcut)
      {
        N_VAddConst(y, ZERO, nhbd_ref);
        generate_nhbd(y, Wrapper.nhbd_fiber, Wrapper.nhbd_partner,
                      Wrapper.mask, params);
      }

      if(t >= scount and scount > 0)
      {
        std::stringstream ss;
        ss << "sindex" << sindex << "tq" << count++;
        save(y, obj, ss.str(), nv_size);
        scount += scount;
      }
    }
    tout = (t + RCONST(8));

    N_VLinearSum(ONE, y, RCONST(-1), y_ref, nv_temp); // Compute y - y_ref

  #ifdef CUDA_TARGET
    // TODO
  #else
    realtype sub_x = NV_Ith_S(nv_temp, sub_index);
    realtype sub_y = NV_Ith_S(nv_temp, sub_index+1);
    NV_Ith_S(nv_temp, sub_index) = 0;
    NV_Ith_S(nv_temp, sub_index + 1) = 0;
  #endif

    equil = N_VMaxNorm(nv_temp);

    realtype adhesion = sqrt(sub_x*sub_x + sub_y*sub_y) / (t - tp);
    realtype term_velocity = sqrt(lambda*lambda + mu*mu);
    bool consider_topsub = params.sub_count == 0 || term_velocity == 0 ?
                              false : true;

    if(consider_topsub and std::abs(term_velocity - adhesion) < movtol)
    {
      if(outcome != 0)
        counter = 0;
      ++counter;
      outcome = 0;
    }else if(equil < movtol and adhesion < movtol)
    {
      if(outcome != 1)
        counter = 0;
      ++counter;
      outcome = 1;
    }else
      counter = 0;

    //std::cout << tp << " " << t << std::endl;
    //std::cout << tout << " " << equil << " " << adhesion << std::endl;
  }
  //clock_t end = clock();
  //std::cout << "Time: " << (end - start)/CLOCKS_PER_SEC << " seconds. "
  //  << "Simulation Time: " << (tout - RCONST(8)) << std::endl;

  if(scount >= 0)
  {
    std::stringstream ss;
    ss << "sindex" << sindex << "tq" << count++;
    save(y, obj, ss.str(), nv_size);
  }

  compute_adhesion(y, params, bot_count, top_count);

  //std::ofstream File("output.json");
  //json::print(File,obj);

  N_VDestroy(y);
  N_VDestroy(y_ref);
  N_VDestroy(nv_temp);
  N_VDestroy(nhbd_ref);
  CVodeFree(&cvode_mem);
  return outcome;
}

int equillibriate_fixed(parameter& params, json::Object& obj, int sindex,
                  realtype scount, uint& bot_count, uint& top_count)
{
  realtype reltol = params.reltol, abstol = params.abstol;
  realtype movtol = params.movtol;
  realtype lambda = params.lambda, mu = params.mu;
  realtype sstep = scount;
  //realtype magnitude = sqrt(lambda*lambda + mu*mu);
  //realtype hmax = magnitude / (params.rmax - params.rcut) / RCONST(1.0001);
  realtype tp = ZERO, t = ZERO, tout = RCONST(8), equil;
  int flag;
  uint count = 0;
  N_Vector y = NULL, nhbd_ref = NULL, y_ref = NULL, nv_temp = NULL;
  N_Vector nv_save = NULL;
  void *cvode_mem = NULL;

  uint particle_count = params.n * params.m + 1;
  uint nv_size = 2*particle_count;
  uint sub_index = 2 * params.n * params.m;
  uint size = params.n * params.m;

  if(params.have_init)
    y = copy_init(params);
  else
    y = generate_init(params);
  y_ref = N_VClone(y);
  nv_temp = N_VClone(y);
  nhbd_ref = N_VClone(y);
  nv_save = N_VClone(y);
  N_VAddConst(y, ZERO, nhbd_ref);

  force_wrapper Wrapper(params, y);

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
  flag = CVodeSetMaxConvFails(cvode_mem, 100);
  flag = CVodeSetMaxStep(cvode_mem, 8);

  //force(RCONST(0), y, y_ref, &Wrapper);
  //N_VPrint_Serial(y_ref);
  //assert(false);


  std::stringstream ss;
  ss << "sindex" << sindex << "tq" << count++;
  save(y, obj, ss.str(), nv_size);

  //clock_t start = clock();
  int counter = 0, outcome = 2;
  realtype ts = 0;
  while(counter < 10)
  {

    N_VAddConst(y, ZERO, y_ref);
    tp = t;
    while(t < tout)
    {
      flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);
      if(flag < 0)
      {
        check_CVode_error(flag);
        outcome = -1;
        N_VPrint_Serial(y);
        exit(-1);
      }

      N_VLinearSum(1, y, -1, nhbd_ref, nv_temp);
      N_VProd(nv_temp, nv_temp, nv_temp); // Square componentwise
      realtype max_dist = ZERO, nmax_dist = ZERO;
      for(int i = 0; i < size; ++i)
      {
        realtype tx = NV_Ith_S(nv_temp, i);
        realtype ty = NV_Ith_S(nv_temp, i + size);
        realtype dist = sqrt(tx*tx + ty*ty);
        if(dist >= max_dist)
        {
          nmax_dist = max_dist;
          max_dist = dist;
        }
      }
      // Also check the upper substrate point.
      realtype tx = NV_Ith_S(nv_temp, 2*size);
      realtype ty = NV_Ith_S(nv_temp, 2*size + 1);
      realtype dist = sqrt(tx*tx + ty*ty);
      if(dist >= max_dist)
      {
        nmax_dist = max_dist;
        max_dist = dist;
      }

      if(max_dist + nmax_dist >= params.rmax - params.rcut)
      {
        N_VAddConst(y, ZERO, nhbd_ref);
        generate_nhbd(y, Wrapper.nhbd_fiber, Wrapper.nhbd_partner,
                      Wrapper.mask, params);
      }

      realtype hu;
      flag = CVodeGetLastStep(cvode_mem, &hu);
      scount = round_nearest_multiple(t - hu, sstep);
      //std::cout << "t - hu: " << (t - hu) << " t: " << t << std::endl;
      //std::cout << "scount: " << scount << std::endl;
      while(scount < t and scount > t - hu)
      {
        printf("%4.4f < %4.4f < %4.4f\n", t-hu, scount, t);
        flag = CVodeGetDky(cvode_mem, scount, 0, nv_save);
        std::stringstream ss1;
        ss1 << "sindex" << sindex << "tq" << count++;
        if(flag >= 0)
          save(nv_save, obj, ss1.str(), nv_size);

        flag = CVodeGetDky(cvode_mem, scount, 1, nv_save);
        std::stringstream ss2;
        ss2 << "sindex" << sindex << "td" << count;
        if(flag >= 0)
          save(nv_save, obj, ss2.str(), nv_size);

        scount += sstep;
      }
    }
    tout = (t + RCONST(8));

    N_VLinearSum(ONE, y, RCONST(-1), y_ref, nv_temp); // Compute y - y_ref

  #ifdef CUDA_TARGET
    // TODO
  #else
    realtype sub_x = NV_Ith_S(nv_temp, sub_index);
    realtype sub_y = NV_Ith_S(nv_temp, sub_index+1);
    NV_Ith_S(nv_temp, sub_index) = 0;
    NV_Ith_S(nv_temp, sub_index + 1) = 0;
  #endif

    equil = N_VMaxNorm(nv_temp);

    realtype adhesion = sqrt(sub_x*sub_x + sub_y*sub_y) / (t - tp);
    realtype term_velocity = sqrt(lambda*lambda + mu*mu);
    bool consider_topsub = params.sub_count == 0 || term_velocity == 0 ?
                              false : true;

    if(consider_topsub and std::abs(term_velocity - adhesion) < movtol)
    {
      if(outcome != 0)
        counter = 0;
      ++counter;
      outcome = 0;
    }else if(equil < movtol and adhesion < movtol)
    {
      if(outcome != 1)
        counter = 0;
      ++counter;
      outcome = 1;
    }else
      counter = 0;

    //std::cout << tp << " " << t << std::endl;
    //std::cout << tout << " " << equil << " " << adhesion << std::endl;
  }
  //clock_t end = clock();
  //std::cout << "Time: " << (end - start)/CLOCKS_PER_SEC << " seconds. "
  //  << "Simulation Time: " << (tout - RCONST(8)) << std::endl;

  if(scount >= 0)
  {
    std::stringstream ss;
    ss << "sindex" << sindex << "tq" << count++;
    save(y, obj, ss.str(), nv_size);
  }

  compute_adhesion(y, params, bot_count, top_count);

  //std::ofstream File("output.json");
  //json::print(File,obj);

  N_VDestroy(y);
  N_VDestroy(y_ref);
  N_VDestroy(nv_temp);
  N_VDestroy(nhbd_ref);
  CVodeFree(&cvode_mem);
  return outcome;
}

realtype round_nearest_multiple(realtype t, realtype multiple)
{
  if(multiple != 0 and t != 0)
  {
    realtype sign = t > ZERO ? ONE : -ONE;
    t *= sign;
    t /= multiple;
    realtype fixed_point = std::ceil(t);
    t = fixed_point*multiple;
    t *= sign;
  }
  return t;
}

void compute_adhesion(N_Vector v, parameter& p, uint& bot_count,uint& top_count)
{
  int size = p.m * p.n;
  for(uint i = 0; i < size; ++i)
  {
    for(uint k = 0; k < p.sub_count; ++k)
    {
      realtype x, xs, y, ys;
      x = NV_Ith_S(v, i);
      y = NV_Ith_S(v, i + size);
      xs = NV_Ith_S(v, 2*size) + k*p.sub_h;
      ys = NV_Ith_S(v, 2*size + 1);

      realtype dist = sqrt((x - xs)*(x - xs) + (y - ys)*(y - ys));
      realtype tol = std::pow(2.0, 1.0/6.0)*p.sigma + 1e-6;
      if(dist <= tol)
      {
        ++top_count;
        break;
      }
    }
  }

  for(uint i = 0; i < size; ++i)
  {
    for(uint k = 0; k < p.osub_count; ++k)
    {
      realtype x, xs, y, ys;
      x = NV_Ith_S(v, i);
      y = NV_Ith_S(v, i + size);
      xs = p.osub + k*p.osub_h;
      ys = ZERO;

      realtype dist = sqrt((x - xs)*(x - xs) + (y - ys)*(y - ys));
      realtype tol = std::pow(2.0, 1.0/6.0)*p.sigma + 1e-6;
      if(dist <= tol)
      {
        ++bot_count;
        break;
      }
    }
  }
}

N_Vector generate_init(parameter& params)
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

N_Vector copy_init(parameter& params)
{
  int n = params.n, m = params.m;
  uint nv_size = TWO*(n*m + ONE);
  N_Vector out = N_VNew_Serial(nv_size);

  for(uint i = 0; i < nv_size; ++i)
  {
    NV_Ith_S(out, i) = params.init[i];
  }

  return out;
}

void handle_kill(int sig)
{
  printf("Caught signal %d\nWriting partial JSON object data.\n",sig);
  (*OBJ_PTR)["grid"] = GRID_PTR;
  std::ofstream File("terminated_output.json");
  json::print(File,*OBJ_PTR);
  exit(1);
}

void check_CVode_error(int flag)
{
  std::cout << "CVode Failed, ";
  switch(flag)
  {
    case CV_MEM_NULL:
    case CV_NO_MALLOC:
    case CV_ILL_INPUT:
      std::cout << "bad memory/malformed input";
      break;
    case CV_ERR_FAILURE:
      std::cout << "too many error test failures during an internal step"
                << " or |h| = h_min";
      break;
    case CV_CONV_FAILURE:
      std::cout << "too many convergence test failures during an internal step"
                << " or |h| = h_min";
      break;
    default:
      std::cout << "TODO: flag = " << flag;
      break;
  }
  std::cout << std::endl;
}
