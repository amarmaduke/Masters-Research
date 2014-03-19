#include "force.h"
#include <stdio.h>

__device__
void position(int& j, int& i, const int ptr, const parameter& p)
{
  j = ptr / p.n;
  i = (ptr % p.n);
}

__global__
void compute_other( value_type* const out_x,
                    value_type* const out_y,
                    value_type* const out_sx,
                    value_type* const out_sy,
                    const value_type* const in_x,
                    const value_type* const in_y,
                    const value_type* const in_s,
                    const parameter p)
{
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  int size = p.n*p.m;

  //if(index >= size)
  //  return;

  int j, i;
  int n = p.n;
  value_type beta = p.beta;
  value_type len = p.len;
  value_type gamma = p.gamma;
  value_type epsilon = p.epsilon;
  value_type sigma = p.sigma;
  value_type* delta = p.delta;

  position(j,i,index,p);

  value_type xp = i == 0 ? delta[j] : in_x[index-1];
  value_type xpp = i == 0 || i == 1 ? delta[j] : in_x[index-2];
  value_type yp = i == 0 ? 0 : in_y[index-1];
  value_type ypp =  i == 0 || i == 1 ?
              ( i == 0 ? -len : 0) : in_y[index-2];

  value_type xn = (index % n) + 1 < n ? in_x[index+1] : NAN;
  value_type xnn = (index % n) + 2 < n ? in_x[index+2] : NAN;
  value_type yn = (index % n) + 1 < n ? in_y[index+1] : NAN;
  value_type ynn = (index % n) + 2 < n ? in_y[index+2] : NAN;

  value_type x = in_x[index];
  value_type y = in_y[index];

  value_type lp = sqrt((xp-xpp)*(xp-xpp) + (yp-ypp)*(yp-ypp));
  value_type l = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp));
  value_type ln = sqrt((xn-x)*(xn-x) + (yn-y)*(yn-y));
  value_type lnn = sqrt((xnn-xn)*(xnn-xn) + (ynn-yn)*(ynn-yn));

  // Bending Force

  value_type xd_f = (xn - x)*(xnn - xn);
  value_type yd_f = (yn - y)*(ynn - yn);
  value_type xd_c = (x - xp)*(xn - x);
  value_type yd_c = (y - yp)*(yn - y);
  value_type xd_b = (xp - xpp)*(x - xp);
  value_type yd_b = (yp - ypp)*(y - yp);

  value_type product_f = xd_f + yd_f;
  value_type product_c = xd_c + yd_c;
  value_type product_b = xd_b + yd_b;

  value_type forward_1 = lnn / (ln * product_f);
  value_type forward_2 = ln*lnn/(product_f*product_f);
  value_type center_1 = ln/(l*product_c);
  value_type center_2 = l/(ln*product_c);
  value_type center_3 = -l*ln/(product_c*product_c);
  value_type backward_1 = lp/(l*product_b);
  value_type backward_2 = lp*l/(product_b*product_b);

  value_type forward_1x = forward_1*(x-xn);
  value_type forward_2x = forward_2*(xnn-xn);
  value_type center_1x = center_1*(x-xp);
  value_type center_2x = center_2*(x-xn);
  value_type center_3x = center_3*(xn-2*x+xp);
  value_type backward_1x = backward_1*(x-xp);
  value_type backward_2x = backward_2*(xpp-xp);

  value_type forward_1y = forward_1*(y-yn);
  value_type forward_2y = forward_2*(ynn-yn);
  value_type center_1y = center_1*(y-yp);
  value_type center_2y = center_2*(y-yn);
  value_type center_3y = center_3*(yn-2*y+yp);
  value_type backward_1y = backward_1*(y-yp);
  value_type backward_2y = backward_2*(ypp-yp);

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

  value_type forward_x = forward_1x + forward_2x;
  value_type center_x = center_1x + center_2x + center_3x;
  value_type backward_x = backward_1x + backward_2x;

  value_type forward_y = forward_1y + forward_2y;
  value_type center_y = center_1y + center_2y + center_3y;
  value_type backward_y = backward_1y + backward_2y;

  value_type bending_x = beta*(forward_x + center_x + backward_x);
  value_type bending_y = beta*(forward_y + center_y + backward_y);

  // Extensible Spring Force

  value_type e_forward = (ln - len)/ln;
  value_type e_backward = (l - len)/l;

  value_type e_forward_x = e_forward*2*(x-xn);
  value_type e_backward_x = e_backward*2*(x-xp);
  value_type e_forward_y = e_forward*2*(y-yn);
  value_type e_backward_y = e_backward*2*(y-yp);

  e_forward_x = isnan(e_forward_x) ? 0 : e_forward_x;
  e_backward_x = isnan(e_backward_x) ? 0 : e_backward_x;
  e_forward_y = isnan(e_forward_y) ? 0 : e_forward_y;
  e_backward_y = isnan(e_backward_y) ? 0 : e_backward_y;

  value_type extensible_x = gamma*(e_forward_x + e_backward_x);
  value_type extensible_y = gamma*(e_forward_y + e_backward_y);

  // Lower substrate vdW pressure

  value_type p1, p2, p4, p5, p7, p8, p11, p13;
  p1 = sigma / y;
  p2 = p1*p1;
  p4 = p2*p2;
  p5 = p4*p1;
  p11 = p5*p5*p1;

  value_type vdW_y = -(PI*epsilon)*(2*p11-4*p5);

  // Upper substrate vdW

  value_type s_x = in_s[0];
  value_type s_y = in_s[1];
  value_type s_vdW_sx = 0, s_vdW_sy = 0, s_vdW_x = 0, s_vdW_y = 0;
  int sub_count = p.sub_count;
  value_type sub_h = p.sub_h;

  for(int k = 0; k < sub_count; ++k)
  {
    value_type x_ = s_x + k*sub_h;
    value_type y_ = s_y;

    value_type xps = x_ - x;
    value_type yps = y_ - y;
    value_type dist = sqrt(xps*xps + yps*yps);

    value_type temp_x = xps/dist;
    value_type temp_y = yps/dist;

    p1 = sigma / dist;
    p2 = p1*p1;
    p4 = p2*p2;
    p7 = p4*p2*p1;
    p8 = p7*p1;
    p13 = p8*p4*p1;
    value_type LJval = -(12*epsilon/sigma)*(p13-p7);

    s_vdW_x = s_vdW_x + LJval*temp_x;
    s_vdW_y = s_vdW_y + LJval*temp_y;

    s_vdW_sx = s_vdW_sx - LJval*temp_x;
    s_vdW_sy = s_vdW_sy - LJval*temp_y;
  }

  // Lower substrate vdW

  value_type os_x = p.osub;
  value_type os_vdW_x = 0, os_vdW_y = 0;
  int osub_count = p.osub_count;
  value_type osub_h = p.osub_h;

  for(int k = 0; k < osub_count; ++k)
  {
    value_type x_ = os_x + k*osub_h;
    value_type y_ = 0;

    value_type xps = x_ - x;
    value_type yps = y_ - y;
    value_type dist = sqrt(xps*xps + yps*yps);

    value_type temp_x = xps/dist;
    value_type temp_y = yps/dist;

    p1 = sigma / dist;
    p2 = p1*p1;
    p4 = p2*p2;
    p7 = p4*p2*p1;
    p8 = p7*p1;
    p13 = p8*p4*p1;
    value_type LJval = -(12*epsilon/sigma)*(p13-p7);

    os_vdW_x = os_vdW_x + LJval*temp_x;
    os_vdW_y = os_vdW_y + LJval*temp_y;
  }

  // Total Force

  value_type total_force_x = -(bending_x + extensible_x) + s_vdW_x + os_vdW_x;
  value_type total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y + os_vdW_y;

  total_force_x = -(bending_x + extensible_x) + s_vdW_x + os_vdW_x;
  total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y + os_vdW_y;

  //printf("index: %d, j: %d, i: %d\nxpp: %f, xp: %f, x: %f, xn: %f, xnn: %f\n ypp: %f, yp: %f, y: %f, yn: %f, ynn: %f\nb_x: %f, e_x: %f, vu_x: %f, vl_x: %f\nb_y: %f, e_y: %f, vp_y: %f, vu_y: %f, vl_y: %f\n",index,j,i,xpp,xp,x,xn,xnn,ypp,yp,y,yn,ynn,bending_x,extensible_x,s_vdW_x,os_vdW_x,bending_y,extensible_y,vdW_y,s_vdW_y,os_vdW_y);

  out_x[index] = total_force_x;
  out_y[index] = total_force_y;
  out_sx[index] = s_vdW_sx;
  out_sy[index] = s_vdW_sy;
}

__device__
value_type2 lennard_jones(value_type2 v, value_type2 v_,
                      int2 idx, int2 idx_, value_type2 acc,
                      value_type sigma, value_type epsilon)
{
  value_type xps = v.x - v_.x;
  value_type yps = v.y - v_.y;

  // Add machine epsilon to prevent 0 / 0 introducing a NaN
  // This is implemented strictly to avoid branching.
  value_type dist = sqrt(xps*xps + yps*yps)
                    + std::numeric_limits<value_type>::epsilon;

  value_type temp_x = xps/dist;
  value_type temp_y = yps/dist;

  value_type p1 = sigma / dist;
  value_type p2 = p1*p1;
  value_type p4 = p2*p2;
  value_type p7 = p4*p2*p1;
  value_type p8 = p7*p1;
  value_type p13 = p8*p4*p1;
  value_type LJval = -(12*epsilon/sigma)*(p13-p7);

  // Condense j == j_ and (i == i_ or i == i_ + 1 or i == i_ - 1) via two
  // always positive continuous functions with zeros only at those points.
  int s1 = ((idx.y - idx_.y)*(idx.y - idx_.y) - 1)*(idx.y - idx_.y)
            *(idx.y - idx_.y)*(idx.y - idx_.y) - 1)*(idx.y - idx_.y);
  int s2 = (idx.x - idx_.x)*(idx.x - idx_.x);

  // Conditional execution instead of branching
  acc.x += ((s1 + s2 == 0)*0) | ((s1 + s2 != 0)*-LJval*temp_x);
  acc.y += ((s1 + s2 == 0)*0) | ((s1 + s2 != 0)*-LJval*temp_y)

  return acc;
}

__device__
value_type2 tile_calculation( value_type2 v, int index, int index_,
                          value_type2 acc, parameter& p)
{
  int i;
  extern __shared__ value_type2 pos[];
  #pragma unroll 4
  for(i = 0; i < blockDim.x; ++i)
  {
    int2 idx, idx_;
    position(idx.x,idx.y,index,p);
    position(idx_.x,idx_.y,index_+i,p);
    acc = lennard_jones(v,pos[i],idx,idx_,acc,p.sigma,p.epsilon);
  }
  return acc;
}

__global__
void compute_n_body(value_type* const out,
                    const value_type* const in,
                    parameter p)
{
  extern __shared__ value_type2 pos[];
  int i, tile, index = blockIdx.x * blockDim.x + threadIdx.x;
  int size = p.n * p.m;
  value_type2 acc = {0.0, 0.0};
  value_type2 v = {in[index], in[index+size]};

  for(i = 0, tile = 0; i < size; i += K, ++tile)
  {
    int idx = tile * blockDim.x + threadIdx.x;
    value_type2 f = {in[idx], in[idx+size]};
    pos[threadIdx.x] = f;
    __syncthreads();
    acc = tile_calculation(v,index,tile*blockDim.x,acc,p);
    __syncthreads();
  }
  out[index] = acc.x;
  out[index+size] = acc.y;
}

void
force_functor::operator() ( const vector_type &x,
                            vector_type &dxdt,
                            const value_type dt)
{
  int size = this->state.n*this->state.m;
  int B = size%K != 0? size/K + 1 : size/K;
  dim3 block_other(B,1,1), thread_other(K,1,1);
  dim3 block_nbody(B,1,1), thread_nbody(K,1,1);

  const value_type* const in = x.data().get();
  value_type* const out = dxdt.data().get();

  cudaStream_t s1, s2;
  cudaEvent_t e1, e2;
  cudaStreamCreate(&s1);
  cudaStreamCreate(&s2);
  cudaEventCreate(&e1);
  cudaEventCreate(&e2);

  thrust::device_ptr<value_type> nbody, substrate;
  nbody = thrust::device_malloc<value_type>(2*size);
  substrate = thrust::device_malloc<value_type>(2*size);
  //cudaMemset(nbody.get(),0,2*size*sizeof(value_type));

  compute_other<<<block_other,thread_other,0>>>
                ( out,out+size,
                  substrate.get(),substrate.get()+size,
                  in,in+size,in+2*size,
                  this->state);
  //cudaEventRecord(e1,s1);

  //cudaDeviceSynchronize();
  compute_n_body<<<block_nbody,thread_nbody,K*sizeof(value_type2)>>>
                  (nbody.get(),in,this->state);

  //cudaDeviceSynchronize();
  //cudaEventSynchronize(e1);

  value_type sub_x = thrust::reduce(substrate,substrate+size);
  value_type sub_y = thrust::reduce(substrate+size,substrate+2*size);

  //cudaDeviceSynchronize();

  //cudaEventSynchronize(e2);
  //cudaEventSynchronize(e1);


  thrust::transform(nbody,nbody+2*size,dxdt.data(),dxdt.data(),thrust::plus<value_type>());

  //cudaDeviceSynchronize();

  dxdt[2*size] = sub_x + this->state.mu;
  dxdt[2*size+1] = sub_y - this->state.lambda;

  cudaStreamDestroy(s1);
  cudaStreamDestroy(s2);
  cudaEventDestroy(e1);
  thrust::device_free(nbody);
  thrust::device_free(substrate);
}
