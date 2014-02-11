#include "force.h"

__global__
void compute_lengths( double * const lens,
                      const double * const in_x,
                      const double * const in_y,
                      const double * const delta,
                      const parameter p)
{
  int j = threadIdx.x;
  int i = threadIdx.y;
  int n = p.n;

  int index = j*n + i;

  double x = in_x[index];
  double y = in_y[index];

  double xp = i == 0 ? delta[j] : in_x[index-1];
  double yp = i == 0 ? 0 : in_y[index-1];

  double l = sqrt((xp-x)*(xp-x)+(yp-y)*(yp-y));

  lens[index] = l;
}

__global__
void compute_fiber_dependent( double * const out_x,
                              double * const out_y,
                              const double * const in_x,
                              const double * const in_y,
                              const double * const lens,
                              const double * const delta,
                              const parameter p)
{
  int j = threadIdx.x; // What fiber we're on
  int i = threadIdx.y; // What particle it is
  int n = p.n;
  double beta = p.beta;
  double len = p.len;
  double gamma = p.gamma;
  double epsilon = p.epsilon;
  double sigma = p.sigma;

  int index = j*n + i;

  double xp = i == 0 ? delta[j] : in_x[index-1];
  double xpp = i == 0 || i == 1 ? delta[j] : in_x[index-2];

  double yp = i == 0 ? 0 : in_y[index-1];
  double ypp =  i == 0 || i == 1 ?
              ( i == 0 ? -len : 0) : in_y[index-2];

  double xn = (index % n) + 1 < n ? in_x[index+1] : NAN;
  double xnn = (index % n) + 2 < n ? in_x[index+2] : NAN;
  double yn = (index % n) + 1 < n ? in_y[index+1] : NAN;
  double ynn = (index % n) + 2 < n ? in_y[index+2] : NAN;

  double l = lens[index];
  double x = in_x[index];
  double y = in_y[index];

  double lp = i == 0 ? sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
              : lens[index-1];
  double ln = (index % n) + 1 < n ? lens[index+1] : NAN;
  double lnn = (index % n) + 2 < n ? lens[index+2] : NAN;

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
  double center_3 = l*ln/(product_c*product_c);
  double backward_1 = lp/(l*product_b);
  double backward_2 = lp*l/(product_b*product_b);

  double forward_1x_t = forward_1*(x-xn);
  double forward_2x_t = forward_2*(xnn-xn);
  double center_1x_t = center_1*(x-xp);
  double center_2x_t = center_2*(x-xn);
  double center_3x_t = center_3*(xn-2*x+xp);
  double backward_1x_t = backward_1*(x-xp);
  double backward_2x_t = backward_2*(xpp-xp);

  double forward_1y_t = forward_1*(y-yn);
  double forward_2y_t = forward_2*(ynn-yn);
  double center_1y_t = center_1*(y-yp);
  double center_2y_t = center_2*(y-yn);
  double center_3y_t = center_3*(yn-2*y+yp);
  double backward_1y_t = backward_1*(y-yp);
  double backward_2y_t = backward_2*(ypp-yp);

  double forward_1x = isnan(forward_1x_t) ? 0 : forward_1x_t;
  double forward_2x = isnan(forward_2x_t) ? 0 : forward_2x_t;
  double center_1x = isnan(center_1x_t) ? 0 : center_1x_t;
  double center_2x = isnan(center_2x_t) ? 0 : center_2x_t;
  double center_3x = isnan(center_3x_t) ? 0 : center_3x_t;
  double backward_1x = isnan(backward_1x_t) ? 0 : backward_1x_t;
  double backward_2x = isnan(backward_2x_t) ? 0 : backward_2x_t;

  double forward_1y = isnan(forward_1y_t) ? 0 : forward_1y_t;
  double forward_2y = isnan(forward_2y_t) ? 0 : forward_2y_t;
  double center_1y = isnan(center_1y_t) ? 0 : center_1y_t;
  double center_2y = isnan(center_2y_t) ? 0 : center_2y_t;
  double center_3y = isnan(center_3y_t) ? 0 : center_3y_t;
  double backward_1y = isnan(backward_1y_t) ? 0 : backward_1y_t;
  double backward_2y = isnan(backward_2y_t) ? 0 : backward_2y_t;

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

  double e_forward_x_t = e_forward*2*(x-xn);
  double e_backward_x_t = e_backward*2*(x-xp);
  double e_forward_y_t = e_forward*2*(y-yn);
  double e_backward_y_t = e_backward*2*(y-yp);

  double e_forward_x = isnan(e_forward_x_t) ? 0 : e_forward_x_t;
  double e_backward_x = isnan(e_backward_x_t) ? 0 : e_backward_x_t;
  double e_forward_y = isnan(e_forward_y_t) ? 0 : e_forward_y_t;
  double e_backward_y = isnan(e_backward_y_t) ? 0 : e_backward_y_t;

  double extensible_x = gamma*(e_forward_x + e_backward_x);
  double extensible_y = gamma*(e_forward_y + e_backward_y);

  // Compute particle-substrate vdW

  double p1 = sigma / y;
  double p2 = p1*p1;
  double p4 = p2*p2;
  double p5 = p4*p1;
  double p11 = p5*p5*p1;

  double vdW_y = -(3.14*epsilon)*(2*p11-4*p5);

  // Compute Total Force

  double total_force_x = -(bending_x + extensible_x);
  double total_force_y = -(bending_y + extensible_y + vdW_y);

  out_x[index] = total_force_x;
  out_y[index] = total_force_y;
}

void compute_n_body_vdw(double * const out_x,
                        double * const out_y,
                        double * const out_s,
                        const double * const in_x,
                        const double * const in_y,
                        const double * const in_s,
                        const parameter p)
{
  int j = threadIdx.x; // What fiber we're on
  int i = threadIdx.y; // What particle it is

  double n = p.n;
  double m = p.m;
  double epsilon = p.epsilon;
  double sigma = p.sigma;
  double sub_h = p.sub_h;
  double sub_count = p.sub_count;

  int index = j*n + i;

  double x = in_x[index];
  double y = in_y[index];

  double vdW_x = 0, vdW_y = 0;
  for(int j_ = 0; j_ < m; ++j_)
  {
    for(int i_ = 0; i_ < n; ++i_)
    {
      if(j == j_ && (i_ == i || i_ + 1 == i || i_ - 1 == i))
        continue;
      int k = j_*n + i_;

      double x_ = in_x[k];
      double y_ = in_y[k];

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

      vdW_x = vdW_x - LJval*temp_x;
      vdW_y = vdW_y - LJval*temp_y;
    }
  }

  double s_x = in_s[0];
  double s_y = in_s[1];
  double s_vdW_x = 0, s_vdW_y = 0;

  for(int k = 0; k < sub_count; ++k)
  {
    double x_ = s_x +k*sub_h;
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

    vdW_x = vdW_x - LJval*temp_x;
    vdW_y = vdW_y - LJval*temp_y;
    s_vdW_x = s_vdW_x + LJval*temp_x;
    s_vdW_y = s_vdW_y + LJval*temp_y;
  }

	atomicAdd(out_s,s_vdW_x);
	atomicAdd(out_s+1,s_vdW_y);

  out_x[index] = vdW_x;
  out_y[index] = vdW_y;
}

__global__
void combine( double * const out_x,
              double * const out_y,
              const double * const f_1x,
              const double * const f_1y,
              const double * const f_2x,
              const double * const f_2y,
              const parameter p)
{
  int j = threadIdx.x; // What fiber we're on
  int i = threadIdx.y; // What particle it is

  double n = p.n;
  int index = j*n + i;

  out_x[index] = f_1x[index] + f_2x[index];
  out_y[index] = f_1y[index] + f_2y[index];
}

thrust::device_ptr<double>
force_functor::operator()(double t, thrust::device_ptr<double> y)
{
  dim3 grid(1,1,1), blocks(this->state.m,this->state.n,1);
  int size = this->state.n * this->state.m;

  thrust::device_ptr<double> out = thrust::device_malloc<double>(2*size+2);
  thrust::device_ptr<double> f_1 = thrust::device_malloc<double>(2*size+2);
  thrust::device_ptr<double> f_2 = thrust::device_malloc<double>(2*size+2);
  thrust::device_ptr<double> lens = thrust::device_malloc<double>(size);

  compute_lengths<<<grid,blocks>>>( lens.get(), y.get(), y.get()+size,
                                    this->state.delta, this->state);
  //std::cout << "functor call: lens:" << std::endl;
	//util::print(lens,size);
	compute_fiber_dependent<<<grid,blocks>>>( f_1.get(), f_1.get()+size, y.get(),
                      y.get()+size, lens.get(), this->state.delta, this->state);
	f_1[2*size] = 0; f_1[2*size+1] = 0;
	f_2[2*size] = 0; f_2[2*size+1] = 0;
	//std::cout << "functor call: fiber_dep:" << std::endl;
	//util::print(f_1,2*size+2);
  compute_n_body_vdw<<<grid,blocks>>>(f_2.get(), f_2.get()+size,
          f_2.get()+2*size, y.get(), y.get()+size, y.get()+2*size, this->state);
	//std::cout << "functor call: n_body:" << std::endl;
	//util::print(f_2,2*size+2);
  thrust::transform(f_1,f_1+2*size+2,f_2,out,thrust::plus<double>());
	
	//std::cout << "functor call: total force:" << std::endl;
	//util::print(out,2*size+2);

  //combine<<<grid,blocks>>>(out_x,out_y,f_1x,f_1y,f_2x,f_2y,p);

  thrust::device_free(f_1);
  thrust::device_free(f_2);
  thrust::device_free(lens);
  return out;
}
