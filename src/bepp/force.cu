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

	double p1, p2, p4, p5, p7, p8, p11, p13;
  p1 = sigma / y;
  p2 = p1*p1;
  p4 = p2*p2;
  p5 = p4*p1;
  p11 = p5*p5*p1;

  double vdW_y = -(PI*epsilon)*(2*p11-4*p5);

  // Compute substrate vdW

  double s_x = in_s[0];
  double s_y = in_s[1];
  double s_vdW_sx = 0, s_vdW_sy = 0, s_vdW_x = 0, s_vdW_y = 0;
  int sub_count = p.sub_count;
  double sub_h = p.sub_h;

  for(int k = 0; k < sub_count; ++k)
  {
    double x_ = s_x + k*sub_h;
    double y_ = s_y;

    double xps = x_ - x;
    double yps = y_ - y;
    double dist = sqrt(xps*xps + yps*yps);

    double temp_x = xps/dist;
    double temp_y = yps/dist;

    p1 = sigma / dist;
    p2 = p1*p1;
    p4 = p2*p2;
    p7 = p4*p2*p1;
    p8 = p7*p1;
    p13 = p8*p4*p1;
    double LJval = -(12*epsilon/sigma)*(p13-p7);

    s_vdW_x = s_vdW_x + LJval*temp_x;
    s_vdW_y = s_vdW_y + LJval*temp_y;

    s_vdW_sx = s_vdW_sx - LJval*temp_x;
    s_vdW_sy = s_vdW_sy - LJval*temp_y;
  }

	// Lower substrate vdW

  double os_x = p.osub;
  double os_vdW_x = 0, os_vdW_y = 0;
  int osub_count = p.osub_count;
  double osub_h = p.osub_h;

  for(int k = 0; k < osub_count; ++k)
  {
    double x_ = os_x + k*osub_h;
    double y_ = 0;

    double xps = x_ - x;
    double yps = y_ - y;
    double dist = sqrt(xps*xps + yps*yps);

    double temp_x = xps/dist;
    double temp_y = yps/dist;

    p1 = sigma / dist;
    p2 = p1*p1;
   	p4 = p2*p2;
    p7 = p4*p2*p1;
    p8 = p7*p1;
    p13 = p8*p4*p1;
    double LJval = -(12*epsilon/sigma)*(p13-p7);

    os_vdW_x = os_vdW_x + LJval*temp_x;
    os_vdW_y = os_vdW_y + LJval*temp_y;
  }

  // Total Force

  double total_force_x = -(bending_x + extensible_x) + s_vdW_x + os_vdW_x;
  double total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y + os_vdW_y;

	total_force_x = -(bending_x + extensible_x) + s_vdW_x + os_vdW_x;
	total_force_y = -(bending_y + extensible_y + vdW_y) + s_vdW_y + os_vdW_y;

	//printf("index: %d, j: %d, i: %d\nxpp: %f, xp: %f, x: %f, xn: %f, xnn: %f\n ypp: %f, yp: %f, y: %f, yn: %f, ynn: %f\nb_x: %f, e_x: %f, vu_x: %f, vl_x: %f\nb_y: %f, e_y: %f, vp_y: %f, vu_y: %f, vl_y: %f\n",index,j,i,xpp,xp,x,xn,xnn,ypp,yp,y,yn,ynn,bending_x,extensible_x,s_vdW_x,os_vdW_x,bending_y,extensible_y,vdW_y,s_vdW_y,os_vdW_y);

	out_x[index] = total_force_x;
	out_y[index] = total_force_y;
	out_sx[index] = s_vdW_sx;
	out_sy[index] = s_vdW_sy;
}

__device__ 
double2 lennard_jones(double2 v, double2 v_, 
											int2 idx, int2 idx_, double2 acc, 
											double sigma, double epsilon)
{
	double xps = v.x - v_.x;
	double yps = v.y - v_.y;
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

	int swtch = (abs(idx.y - idx_.y) - 1)*(idx.y - idx_.y);

	acc.x += idx.x == idx_.x and swtch == 0? 0 : -LJval*temp_x;
	acc.y += idx.x == idx_.x and swtch == 0? 0 : -LJval*temp_y;

	//printf("block: %d, thread: %d, j: %d, j_: %d, i: %d, i_: %d\nx: %f, y: %f, x_: %f, y_: %f\ndist: %f, acc.x: %f, acc.y: %f\n",blockIdx.x,threadIdx.x,idx.x,idx_.x,idx.y,idx_.y,v.x,v.y,v_.x,v_.y,dist,acc.x,acc.y);
	//assert(not isnan(acc.x) and not isnan(acc.y));

	return acc;
}

__device__ 
double2 tile_calculation(	double2 v, int index, int index_,
													double2 acc, parameter& p)
{
	int i;
	extern __shared__ double2 pos[];
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
void compute_n_body2(double* const out, const double* const in, parameter p)
{
	extern __shared__ double2 pos[];
	int i, tile, index = blockIdx.x * blockDim.x + threadIdx.x;
	int size = p.n * p.m;
	double2 acc = {0.0, 0.0};
	double2 v = {in[index], in[index+size]};

	for(i = 0, tile = 0; i < size; i += K, ++tile)
	{
		int idx = tile * blockDim.x + threadIdx.x;
		double2 f = {in[idx], in[idx+size]};
		pos[threadIdx.x] = f;
		__syncthreads();
		acc = tile_calculation(v,index,tile*blockDim.x,acc,p);
		__syncthreads();
	}
	out[index] = acc.x;
	out[index+size] = acc.y;
}

__global__
void compute_n_body(double* const out_x,
                    double* const out_y,
                    const double* const in_x,
                    const double* const in_y,
                    const parameter p)
{
	//__shared__ double pos_x[K];
	//__shared__ double pos_y[K];
  int row = blockDim.x * blockIdx.x + threadIdx.x;
	int size = p.n*p.m;
	double sigma = p.sigma;
  double epsilon = p.epsilon;

  int j, i, j_, i_;
	double x, x_, y, y_;

	for(int b = 0; b < gridDim.x; ++b)
	{
		int column = blockDim.x * b;
		//__syncthreads();
		/*
		if(column+threadIdx.x < size)
		{
			pos_x[threadIdx.x] = in_x[column+threadIdx.x];
			pos_y[threadIdx.x] = in_y[column+threadIdx.x];
			//printf("pos_x: %f, pos_y: %f, index: %d\n",pos_x[threadIdx.x],pos_y[threadIdx.x],threadIdx.x);
			//assert(pos_y[threadIdx.x] < 100);
		}
		*/
		//__syncthreads();
		
		//#pragma unroll
		for(int k = 0; k < blockDim.x; ++k)
		{
 	 		position(j,i,row,p);
 	 		position(j_,i_,column+k,p);
		
			if(row < size and column+k < size)
			{
			x = in_x[row];
			y = in_y[row];
			x_ = in_x[column+k];
			y_ = in_y[column+k];
			
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

			vdW_x = (j == j_ and (i == i_ or i+1 == i_ or i-1 == i_))? 0 : vdW_x;
			vdW_y = (j == j_ and (i == i_ or i+1 == i_ or i-1 == i_))? 0 : vdW_y;

			printf("x: %f, y: %f, x_: %f, y_: %f\nrow: %d, col: %d\nj: %d, i: %d, j_: %d, i_: %d\nvdW_x: %f, vdW_y: %f\nsize: %d\n\n",x,y,x_,y_,row,column+k,j,i,j_,i_,vdW_x,vdW_y,size);
			//assert(x < 100 and y < 100);
			assert(not isnan(vdW_x) and not isnan(vdW_y));
			out_x[row] += vdW_x;
			out_y[row] += vdW_y;
			}
		}
		__syncthreads();
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
	cudaMemset(nbody.get(),0,2*size*sizeof(value_type));

	/*
	std::cout << "In: " << std::endl;
	double* t1 = new double[2*size+2];
	cudaMemcpy(t1,in,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size+2; ++i)
	{
		printf("i: %d, %4.11f\n",i,t1[i]);
	}
	printf("\n");
	*/

	compute_other<<<block_other,thread_other,0>>>
                (out,out+size,substrate.get(),substrate.get()+size,in,in+size,in+2*size,this->state);
	//cudaEventRecord(e1,s1);

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

	//compute_n_body<<<block_nbody,thread_nbody,0,s1>>>
  //              (nbody,nbody+size,in,in+size,this->state);
	
	//cudaDeviceSynchronize();
	compute_n_body2<<<block_nbody,thread_nbody,K*sizeof(double2)>>>
									(nbody.get(),in,this->state);
	//cudaEventRecord(e2,s2);
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
	
	//cudaDeviceSynchronize();
	//cudaEventSynchronize(e1);

	double sub_x = thrust::reduce(substrate,substrate+size);
	double sub_y = thrust::reduce(substrate+size,substrate+2*size);
	
	//cudaDeviceSynchronize();

	/*
	double* t1 = new double[2*size];
	double* t2 = new double[2*size];
	double* t = new double[2*size+2];
	cudaMemcpy(t1,nbody.get(),sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	cudaMemcpy(t2,substrate.get(),sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	cudaMemcpy(t,out,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size; ++i)
	{
		printf("i: %d, nbody: %f, sub: %f, out: %f\n",i,t1[i],t2[i],t[i]);
		assert(not isnan(t1[i]) and not isnan(t2[i]) and not isnan(t[i]));
	}
	*/

	/*
	std::cout << "Start:" << std::endl;
	double* t = new double[2*size];
	cudaMemcpy(t,nbody,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	double* t1 = new double[2*size];
	cudaMemcpy(t1,out,sizeof(double)*2*size,cudaMemcpyDeviceToHost);
	std::cout << "[";
	for(int i = 0; i < 2*size; ++i)
	{
		std::cout << t1[i] << ",";
	}
	std::cout << "]" << std::endl;
	*/

	//cudaEventSynchronize(e2);
	//cudaEventSynchronize(e1);

	//combine<<<1,2*size>>>(out,nbody);
	thrust::transform(nbody,nbody+2*size,dxdt.data(),dxdt.data(),thrust::plus<double>());

	//cudaDeviceSynchronize();

  dxdt[2*size] = sub_x + this->state.mu;
  dxdt[2*size+1] = sub_y - this->state.lambda;

/*
	std::cout << "In: " << std::endl;
	double* t1 = new double[2*size+2];
	cudaMemcpy(t1,in,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size+2; ++i)
	{
		printf("i: %d, %4.11f\n",i,t1[i]);
		assert(not isnan(t1[i]));
	}
	printf("\n");

	std::cout << "Force:" << std::endl;
	double* t2 = new double[2*size+2];
	cudaMemcpy(t2,out,sizeof(double)*(2*size+2),cudaMemcpyDeviceToHost);
	for(int i = 0; i < 2*size+2; ++i)
	{
		printf("i: %d, %4.11f\n",i,t2[i]);
		assert(not isnan(t1[i]));
	}
	printf("\n");
*/
	//assert(false);
	//assert(not isnan(t2[0]));
  
	cudaStreamDestroy(s1);
 	cudaStreamDestroy(s2);
	cudaEventDestroy(e1);
	thrust::device_free(nbody);
	thrust::device_free(substrate);
}
