#include "force.h"
#include <stdio.h>

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
                    double* const out_sx,
										double* const out_sy,
                    const double* const in_x,
                    const double* const in_y,
                    const double* const in_s,
                    const parameter p)
{
  int index = blockDim.x * blockIdx.x + threadIdx.x;
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

	//printf("xpp: %f, xp: %f, x: %f, xn: %f, xnn: %f\n ypp: %f, yp: %f, y: %f, yn: %f, ynn: %f\n lp: %f, l: %f, ln: %f, lnn: %f\n",xpp,xp,x,xn,xnn,ypp,yp,y,yn,ynn,lp,l,ln,lnn);

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

    double xps = x_ - x;
    double yps = y_ - y;
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

    s_vdW_x = s_vdW_x + LJval*temp_x;
    s_vdW_y = s_vdW_y + LJval*temp_y;

    s_vdW_sx = s_vdW_sx - LJval*temp_x;
    s_vdW_sy = s_vdW_sy - LJval*temp_y;
  }

  // Total Force

  double total_force_x = -(bending_x + extensible_x) + s_vdW_x;
  double total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y;

	//printf("b_x: %f, e_x: %f, e_v: %f\nb_y: %f, e_y: %f, e_v: %f, e_vs: %f\n",bending_x,extensible_x,s_vdW_x,bending_y,extensible_y,s_vdW_y,vdW_y);

	out_x[index] = total_force_x;
	out_y[index] = total_force_y;
	out_sx[index] = s_vdW_sx;
	out_sy[index] = s_vdW_sy;
}

__global__
void compute_n_body(double* const out_x,
                    double* const out_y,
                    const double* const in_x,
                    const double* const in_y,
                    const parameter p)
{
	__shared__ double pos_x[K];
	__shared__ double pos_y[K];
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
	const int size = p.n*p.m;
	const double sigma = p.sigma;
  const double epsilon = p.epsilon;

  int j, i, j_, i_;
	double x, x_, y, y_;
	
	if(row < size)
	{
		x = in_x[row];
		y = in_y[row];
	}

	for(int b = 0; b < blockDim.x; ++b)
	{
		const int column = blockDim.x * b;
		if(column+threadIdx.x < size)
		{
			pos_x[threadIdx.x] = in_x[column+threadIdx.x];
			pos_y[threadIdx.x] = in_y[column+threadIdx.x];
		}
		__syncthreads();
		
		#pragma unroll
		for(int k = 0; k < K; ++k)
		{
 	 		position(j,i,row,p);
 	 		position(j_,i_,column+k,p);
		
			if(row < size and column+k < size and (j != j_ or (i_ != i and i_ + 1 != i and i_ - 1 != i)))
    	{
				x_ = pos_x[k];
				y_ = pos_y[k];
				
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

				//printf("x: %f, y: %f, x_: %f, y_: %f\nrow: %d, col: %d\nj: %d, i: %d, j_: %d, i_: %d\nvdW_x: %f, vdW_y: %f\nsize: %d\n\n",x,y,x_,y_,row,column+k,j,i,j_,i_,vdW_x,vdW_y,size);
				//assert(x < 100);
				//assert(not isnan(vdW_x));
				out_x[row] += vdW_x;
				out_y[row] += vdW_y;
  		}
		}
	}
}

__global__
void combine(double* const out, const double* const in)
{
	int index = threadIdx.x;
	out[index] += in[index];
}

void
force_functor::operator() ( const vector_type &x,
                            vector_type &dxdt,
                            const value_type dt)
{
  int size = this->state.n*this->state.m;
	int total_size = 2*size + 2;
	int B = size%K != 0? size/K + 1 : size/K;
  dim3 block_other(B,1,1), thread_other(K,1,1);
  dim3 block_nbody(B,1,1), thread_nbody(K,1,1);

  const value_type* const in = x.data().get();
  value_type* const out = dxdt.data().get();

  cudaStream_t s1, s2;
	cudaEvent_t e1;
	cudaStreamCreate(&s1);
	cudaStreamCreate(&s2);
	cudaEventCreate(&e1);

  value_type* nbody, *substrate;
	cudaMalloc(&nbody,2*size*sizeof(value_type));
	cudaMalloc(&substrate,2*size*sizeof(value_type));
	cudaMemset(nbody,0,2*size*sizeof(value_type));

	compute_other<<<block_other,thread_other,0,s1>>>
                (out,out+size,substrate,substrate+size,in,in+size,in+2*size,this->state);
	cudaEventRecord(e1,s1);

		/*
		std::cout << "In: " << std::endl;
		double* t2 = new double[2*size+2];
		cudaMemcpy(t2,in,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
		for(int j = 0; j < 2*size+2; ++j)
		{
			std::cout << t2[j] << " ";
		}
		std::cout << std::endl << std::endl;
		*/

	compute_n_body<<<block_nbody,thread_nbody,0,s2>>>
                (nbody,nbody+size,in,in+size,this->state);
		
		
		/*
		std::cout << "Ptrs: K*i: " << (K*i) << std::endl;
		std::cout << "outptr: " << out_ptr << " inptr: " << in_ptr << std::endl;
		std::cout << "nbody: " << nbody << " in: " << in << std::endl << std::endl;
		std::cout << "In: " << std::endl;
		double* t2 = new double[2*size+2];
		cudaMemcpy(t2,in,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
		for(int j = 0; j < 2*size+2; ++j)
		{
			std::cout << t2[j] << " ";
		}
		std::cout << std::endl << std::endl;

		std::cout << "Nbody: " << std::endl;
		double* t1 = new double[2*size];
		cudaMemcpy(t1,nbody,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
		for(int j = 0; j < 2*size; ++j)
		{
			std::cout << t1[j] << " ";
		}
		std::cout << std::endl << std::endl;
	
		free(t1);
		free(t2);
		*/
	
	cudaEventSynchronize(e1);

	double sub_x = thrust::reduce(thrust::device,substrate,substrate+size);
	double sub_y = thrust::reduce(thrust::device,substrate+size,substrate+2*size);
	
	cudaDeviceSynchronize();

	/*
	std::cout << "Start:" << std::endl;
	double* t = new double[2*size];
	cudaMemcpy(t,nbody,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	double* t1 = new double[2*size];
	cudaMemcpy(t1,out,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size; ++i)
	{
		std::cout << t[i] << " " << t1[i] << std::endl;
	}*/

	combine<<<1,2*size>>>(out,nbody);
	//thrust::transform(nbody,nbody+2*size,dxdt.begin(),dxdt.begin(),thrust::plus<double>());

/*
	std::cout << "Next:" << std::endl;
	double* t2 = new double[2*size];
	cudaMemcpy(t2,out,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size; ++i)
	{
		std::cout << t2[i] << std::endl;
	}*/

  dxdt[2*size] = sub_x + this->state.mu;
  dxdt[2*size+1] = sub_y - this->state.lambda;
  
	cudaStreamDestroy(s1);
 	cudaStreamDestroy(s2);
	cudaEventDestroy(e1);
	cudaFree(nbody);
	cudaFree(substrate);
}
