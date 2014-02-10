#include <iostream>
#include <vector>
#include <stdio.h>

#include "force.h"
#include "../solver/dorpi.h"

int main()
{
  double t = 0;
  double step = .0001;
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
  o.tolerance = .00000001;

  time_t start = time(0);
  thrust::host_vector<thrust::device_ptr<double> > grid
      = dorpi::solve(F, t, init, 2*size+2, o);
  time_t end = time(0);
  double time = difftime(end,start);

  for(int i = 0; i < grid.size(); ++i)
  {
    for(int j = 0; j < 2*size+2; ++j)
    {
      std::cout << grid[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "time: " << time << std::endl;

  thrust::device_free(init);
  thrust::device_free(d_delta);
  return 0;
}
