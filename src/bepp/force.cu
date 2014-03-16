#include "force.h"

__device__
double atomicAdd(double * const address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(val +
                      __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

__device__
void position(int& j, int& i, const int ptr, const parameter& p)
{
  j = ptr / p.n;
  i = (ptr % p.n);
}

__global__
void compute_other( double* const out_x,
                    double* const out_y,
                    double* const out_s,
                    const double* const in_x,
                    const double* const in_y,
                    const double* const in_s,
                    const parameter p)
{
  int index = blockIdx.x + threadIdx.x;
  int size = p.n*p.m;
  if(index >= size)
    return;

  int j, i;
  int n = p.n;
  double beta = p.beta;
  double len = p.len;
  double gamma = p.gamma;
  double epsilon = p.epsilon;
  double sigma = p.sigma;
  double* delta = p.delta;

  position(j,i,index,p);

  double xp = i == 0 ? delta[j] : in_x[index-1];
  double xpp = i == 0 || i == 1 ? delta[j] : in_x[index-2];
  double yp = i == 0 ? 0 : in_y[index-1];
  double ypp =  i == 0 || i == 1 ?
              ( i == 0 ? -len : 0) : in_y[index-2];

  double xn = (index % n) + 1 < n ? in_x[index+1] : NAN;
  double xnn = (index % n) + 2 < n ? in_x[index+2] : NAN;
  double yn = (index % n) + 1 < n ? in_y[index+1] : NAN;
  double ynn = (index % n) + 2 < n ? in_y[index+2] : NAN;

  double x = in_x[index];
  double y = in_y[index];

  double lp = sqrt((xp-xpp)*(xp-xpp) + (yp-ypp)*(yp-ypp));
  double l = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp));
  double ln = sqrt((xn-x)*(xn-x) + (yn-y)*(yn-y));
  double lnn = sqrt((xnn-xn)*(xnn-xn) + (ynn-yn)*(ynn-yn));

  // Compute Bending Force

  double xd_f = (xn - x)*(xnn - xn);
  double yd_f = (yn - y)*(ynn - yn);
  double xd_c = (x - xp)*(xn - x);
  double yd_c = (y - yp)*(yn - y);
  double xd_b = (xp - xpp)*(x - xp);
  double yd_b = (yp - ypp)*(y - yp);

  double product_f = xd_f + yd_f;
  double product_c = xd_c + yd_c;
  double product_b = xd_b + yd_b;

  double forward_1 = lnn / (ln * product_f);
  double forward_2 = ln*lnn/(product_f*product_f);
  double center_1 = ln/(l*product_c);
  double center_2 = l/(ln*product_c);
  double center_3 = -l*ln/(product_c*product_c);
  double backward_1 = lp/(l*product_b);
  double backward_2 = lp*l/(product_b*product_b);

  double forward_1x = forward_1*(x-xn);
  double forward_2x = forward_2*(xnn-xn);
  double center_1x = center_1*(x-xp);
  double center_2x = center_2*(x-xn);
  double center_3x = center_3*(xn-2*x+xp);
  double backward_1x = backward_1*(x-xp);
  double backward_2x = backward_2*(xpp-xp);

  double forward_1y = forward_1*(y-yn);
  double forward_2y = forward_2*(ynn-yn);
  double center_1y = center_1*(y-yp);
  double center_2y = center_2*(y-yn);
  double center_3y = center_3*(yn-2*y+yp);
  double backward_1y = backward_1*(y-yp);
  double backward_2y = backward_2*(ypp-yp);

  forward_1x = isnan(forward_1x) ? 0 : forward_1x;
  forward_2x = isnan(forward_2x) ? 0 : forward_2x;
  center_1x = isnan(center_1x) ? 0 : center_1x;
  center_2x = isnan(center_2x) ? 0 : center_2x;
  center_3x = isnan(center_3x) ? 0 : center_3x;
  backward_1x = isnan(backward_1x) ? 0 : backward_1x;
  backward_2x = isnan(backward_2x) ? 0 : backward_2x;

  forward_1y = isnan(forward_1y) ? 0 : forward_1y;
  forward_2y = isnan(forward_2y) ? 0 : forward_2y;
  center_1y = isnan(center_1y) ? 0 : center_1y;
  center_2y = isnan(center_2y) ? 0 : center_2y;
  center_3y = isnan(center_3y) ? 0 : center_3y;
  backward_1y = isnan(backward_1y) ? 0 : backward_1y;
  backward_2y = isnan(backward_2y) ? 0 : backward_2y;

  double forward_x = forward_1x + forward_2x;
  double center_x = center_1x + center_2x + center_3x;
  double backward_x = backward_1x + backward_2x;

  double forward_y = forward_1y + forward_2y;
  double center_y = center_1y + center_2y + center_3y;
  double backward_y = backward_1y + backward_2y;

  double bending_x = beta*(forward_x + center_x + backward_x);
  double bending_y = beta*(forward_y + center_y + backward_y);

  // Compute Extensible Spring Force

  double e_forward = (ln - len)/ln;
  double e_backward = (l - len)/l;

  double e_forward_x = e_forward*2*(x-xn);
  double e_backward_x = e_backward*2*(x-xp);
  double e_forward_y = e_forward*2*(y-yn);
  double e_backward_y = e_backward*2*(y-yp);

  e_forward_x = isnan(e_forward_x) ? 0 : e_forward_x;
  e_backward_x = isnan(e_backward_x) ? 0 : e_backward_x;
  e_forward_y = isnan(e_forward_y) ? 0 : e_forward_y;
  e_backward_y = isnan(e_backward_y) ? 0 : e_backward_y;

  double extensible_x = gamma*(e_forward_x + e_backward_x);
  double extensible_y = gamma*(e_forward_y + e_backward_y);

  // Compute particle-substrate vdW

  double p1 = sigma / y;
  double p2 = p1*p1;
  double p4 = p2*p2;
  double p5 = p4*p1;
  double p11 = p5*p5*p1;

  double vdW_y = -(PI*epsilon)*(2*p11-4*p5);

  // Compute substrate vdW

  double s_x = in_s[0];
  double s_y = in_s[1];
  double s_vdW_sx = 0, s_vdW_sy = 0, s_vdW_x = 0, s_vdW_y = 0;
  int sub_count = p.sub_count;
  int sub_h = p.sub_h;

  for(int k = 0; k < sub_count; ++k)
  {
    double x_ = s_x + k*sub_h;
    double y_ = s_y;

    double xps = x - x_;
    double yps = y - y_;
    double dist = sqrt(xps*xps + yps*yps);

    double temp_x = xps/dist;
    double temp_y = yps/dist;

    double p1 = sigma / dist;
    double p2 = p1*p1;
    double p4 = p2*p2;
    double p7 = p4*p2*p1;
    double p8 = p7*p1;
    double p13 = p8*p4*p1;
    double LJval = -(12*epsilon/sigma)*(p13-p7);

    s_vdW_x = s_vdW_x - LJval*temp_x;
    s_vdW_y = s_vdW_y - LJval*temp_y;

    s_vdW_sx = s_vdW_sx + LJval*temp_x;
    s_vdW_sy = s_vdW_sy + LJval*temp_y;
  }

  // Total Force

  double total_force_x = -(bending_x + extensible_x) + s_vdW_x;
  double total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y;

  atomicAdd(out_x+index,total_force_x);
  atomicAdd(out_y+index,total_force_y);
  atomicAdd(out_s,s_vdW_sx);
  atomicAdd(out_s+1,s_vdW_sy);
}

__global__
void compute_n_body(double* const out_x,
                    double* const out_y,
                    const double* const in_x,
                    const double* const in_y,
                    const parameter p)
{
  const int row = blockIdx.y + threadIdx.y;
  const int col = blockIdx.x + threadIdx.x;
  const int size = p.n*p.m;
  const double sigma = p.sigma;
  const double epsilon = p.epsilon;

  int j, i, j_, i_;
  position(j,i,row,p);
  position(j_,i_,col,p);

  if(row >= size or col >= size
      or (j == j_ and (i_ == i or i_ + 1 == i or i_ - 1 == i)))
  {
	}else
  {
    double x = in_x[row];
    double y = in_y[row];
    double x_ = in_x[col];
    double y_ = in_y[col];

    double xps = x - x_;
    double yps = y - y_;
    double dist = sqrt(xps*xps + yps*yps);

    double temp_x = xps/dist;
    double temp_y = yps/dist;

    double p1 = sigma / dist;
    double p2 = p1*p1;
    double p4 = p2*p2;
    double p7 = p4*p2*p1;
    double p8 = p7*p1;
    double p13 = p8*p4*p1;
    double LJval = -(12*epsilon/sigma)*(p13-p7);

    double vdW_x = -LJval*temp_x;
    double vdW_y = -LJval*temp_y;

    atomicAdd(out_x+row,vdW_x);
    atomicAdd(out_y+row,vdW_y);
  }
}

void
force_functor::operator() ( const vector_type &x,
                            vector_type &dxdt,
                            const value_type dt)
{
  int size = this->state.n * this->state.m;
  int B = (size+1)/K;
  dim3 block_other(B,1,1), thread_other(K,1,1);
  dim3 block_nbody(B,B,1), thread_nbody(K,K,1);

  const double* in = x.data().get();
  double* out = dxdt.data().get();

  cudaStream_t s1, s2;
  cudaStreamCreate(&s1);
  cudaStreamCreate(&s2);

  compute_other<<<block_other,thread_other,0,s1>>>
                (out,out+size,out+2*size,in,in+size,in+2*size,this->state);
  compute_n_body<<<block_nbody,thread_nbody,0,s2>>>
                (out,out+size,in,in+size,this->state);

  dxdt[2*size] += this->state.mu;
  dxdt[2*size+1] -= this->state.lambda;

  cudaStreamDestroy(s1);
  cudaStreamDestroy(s2);
}
