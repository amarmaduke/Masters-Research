
#include "force.h"

extern "C" int force(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  N_VScale(ZERO, udot, udot);

  force_wrapper& wp = *((force_wrapper*)(user_data));
  parameter params = wp.params;
  int n = params.n, m = params.m, size = n*m;

  N_Vector vtemp = N_VClone(u);
  N_VLinearSum(1, u, -1, wp.u_ref, vtemp);
  N_VProd(vtemp, vtemp, vtemp); // Square componentwise
  realtype max_dist = ZERO;
  for(int i = 0; i < size; ++i)
  {
    realtype x = NV_Ith_S(vtemp, i);
    realtype y = NV_Ith_S(vtemp, i + size);
    realtype dist = sqrt(x*x + y*y);
    max_dist = std::max(max_dist, dist);
  }
  // Also check the upper substrate point.
  realtype sx = NV_Ith_S(vtemp, 2*size);
  realtype sy = NV_Ith_S(vtemp, 2*size + 1);
  realtype sdist = sqrt(sx*sx + sy*sy);
  max_dist = std::max(max_dist, sdist);
  N_VDestroy(vtemp);

  if(max_dist >= params.radius/2)
  {
    N_VAddConst(u, 0, wp.u_ref);
    generate_nhbd(u, wp.nhbd_fiber, wp.nhbd_partner, wp.mask, params);
  }

  // Moving substrate load
  NV_Ith_S(udot, 2*size) += params.mu;
  NV_Ith_S(udot, 2*size+1) -= params.lambda;

  // van der Waals fiber-fiber, fiber-substrate
  int f_idx, p_idx, p_xindex, p_yindex;
  bool is_lower_substrate = false;
  realtype xps, yps, dist, temp_x, temp_y;
  realtype p1, p2, p4, p7, p8, p13, LJval;
  realtype epsi = params.epsilon;
  for(int i = 0, b = 0; i < wp.nhbd_fiber.size(); ++i, b+=2)
  {
    is_lower_substrate = false;
    realtype x_f, y_f, x_p, y_p;
    f_idx = wp.nhbd_fiber[i];
    p_idx = wp.nhbd_partner[i];

    x_f = NV_Ith_S(u, f_idx);
    y_f = NV_Ith_S(u, f_idx + size);

    if(wp.mask[b] and wp.mask[b+1]) // Fiber-Fiber
    {
      x_p = NV_Ith_S(u, p_idx);
      y_p = NV_Ith_S(u, p_idx + size);
      epsi = params.epsilon;
      p_xindex = p_idx;
      p_yindex = p_xindex + size;
    }else if(wp.mask[b] and not wp.mask[b+1]) // Fiber-Upper
    {
      x_p = NV_Ith_S(u, 2*size) + p_idx*params.sub_h;
      y_p = NV_Ith_S(u, 2*size + 1);
      epsi = params.epsilon_top;
      p_xindex = 2*size;
      p_yindex = p_xindex + 1;
    }else if(not wp.mask[b] and wp.mask[b+1]) // Fiber-Lower
    {
      x_p = params.osub + p_idx*params.osub_h;
      y_p = ZERO;
      epsi = params.epsilon_bottom;
      is_lower_substrate = true;
    }

    xps = x_f - x_p;
    yps = y_f - y_p;
    dist = sqrt(xps*xps + yps*yps);
    temp_x = xps/dist;
    temp_y = yps/dist;

    p1 = params.sigma/dist;
    p2 = p1*p1; p4 = p2*p2; p8 = p4*p4;
    p7 = p4*p2*p1; p13 = p8*p4*p1;
    LJval = -(RCONST(12)*epsi/params.sigma)*(p13 - p7);

    NV_Ith_S(udot, f_idx) -= LJval*temp_x;
    NV_Ith_S(udot, f_idx + size) -= LJval*temp_y;
    if(not is_lower_substrate)
    {
      NV_Ith_S(udot, p_xindex) += LJval*temp_x;
      NV_Ith_S(udot, p_yindex) += LJval*temp_y;
    }
  }

  realtype total_fiber_force_x, total_fiber_force_y;
  realtype xpp, xp, x, xn, xnn, ypp, yp, y, yn, ynn, lp, l, ln, lnn;

  for(int j = 0; j < m; ++j)
  {
    for(int i = 0; i < n; ++i)
    {
      // Setup required points
      int xindex = i + j*n, yindex = xindex + size;
      xpp = i == 0 or i == 1 ? params.delta[j] : NV_Ith_S(u, xindex-2);
      xp = i == 0 ? params.delta[j] : NV_Ith_S(u, xindex-1);
      ypp = i == 0 ? -params.len : (i == 1 ? ZERO : NV_Ith_S(u, yindex-2));
      yp = i == 0 ? ZERO : NV_Ith_S(u, yindex-1);
      xn = (xindex % n) + 1 < n ? NV_Ith_S(u, xindex+1) : NAN;
      xnn = (xindex % n) + 2 < n ? NV_Ith_S(u, xindex+2) : NAN;
      yn = (xindex % n) + 1 < n ? NV_Ith_S(u, yindex+1) : NAN;
      ynn = (xindex % n) + 2 < n ? NV_Ith_S(u, yindex+2) : NAN;
      x = NV_Ith_S(u, xindex);
      y = NV_Ith_S(u, yindex);

      // Compute required lengths
      lp = sqrt((xp-xpp)*(xp-xpp) + (yp-ypp)*(yp-ypp));
      l = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp));
      ln = sqrt((xn-x)*(xn-x) + (yn-y)*(yn-y));
      lnn = sqrt((xnn-xn)*(xnn-xn) + (ynn-yn)*(ynn-yn));

      // Compute Bending
      realtype xd_f, yd_f, xd_c, yd_c, xd_b, yd_b;
      xd_f = (xn - x)*(xnn - xn);
      yd_f = (yn - y)*(ynn - yn);
      xd_c = (x - xp)*(xn - x);
      yd_c = (y - yp)*(yn - y);
      xd_b = (xp - xpp)*(x - xp);
      yd_b = (yp - ypp)*(y - yp);

      realtype product_f, product_c, product_b;
      product_f = xd_f + yd_f;
      product_c = xd_c + yd_c;
      product_b = xd_b + yd_b;

      realtype b_3_t1, b_3_t2, b_3_b, b_3x, b_3y;
      b_3_t1 = FOUR*(lnn/ln)*product_f;
      b_3_t2 = FOUR*(lnn*ln);
      b_3_b = (lnn*ln + product_f)*(lnn*ln+product_f);
      b_3x = (b_3_t1*(x-xn) - b_3_t2*(xn-xnn))/b_3_b;
      b_3y = (b_3_t1*(y-yn) - b_3_t2*(yn-ynn))/b_3_b;

      realtype b_2_t1x, b_2_t2x, b_2_t1y, b_2_t2y, b_2_b, b_2x, b_2y;
      b_2_t1x = FOUR*(l/ln*(x-xn) + ln/l*(x-xp))*product_c;
      b_2_t2x = FOUR*ln*l*(xp-TWO*x+xn);
      b_2_t1y = FOUR*(l/ln*(y-yn) + ln/l*(y-yp))*product_c;
      b_2_t2y = FOUR*ln*l*(yp-TWO*y+yn);
      b_2_b = (ln*l + product_c)*(ln*l + product_c);
      b_2x = (b_2_t1x - b_2_t2x)/b_2_b;
      b_2y = (b_2_t1y - b_2_t2y)/b_2_b;

      realtype b_1_t1, b_1_t2, b_1_b, b_1x, b_1y;
      b_1_t1 = FOUR*(lp/l)*product_b;
      b_1_t2 = FOUR*l*lp;
      b_1_b = (l*lp + product_b)*(l*lp + product_b);
      b_1x = (b_1_t1*(x-xp) - b_1_t2*(xp-xpp))/b_1_b;
      b_1y = (b_1_t1*(y-yp) - b_1_t2*(yp-ypp))/b_1_b;

      b_1x = isnan(b_1x) ? ZERO : b_1x;
      b_2x = isnan(b_2x) ? ZERO : b_2x;
      b_3x = isnan(b_3x) ? ZERO : b_3x;

      b_1y = isnan(b_1y) ? ZERO : b_1y;
      b_2y = isnan(b_2y) ? ZERO : b_2y;
      b_3y = isnan(b_3y) ? ZERO : b_3y;

      realtype bending_x = params.beta*(b_1x + b_2x + b_3x);
      realtype bending_y = params.beta*(b_1y + b_2y + b_3y);

      // Compute Extensible Spring
      realtype e_forward, e_backward;
      e_forward = (ln - params.len)/ln;
      e_backward = (l - params.len)/l;

      realtype e_forwardx, e_forwardy, e_backwardx, e_backwardy;
      e_forwardx = e_forward*TWO*(x-xn);
      e_backwardx = e_backward*TWO*(x-xp);
      e_forwardy = e_forward*TWO*(y-yn);
      e_backwardy = e_backward*TWO*(y-yp);

      e_forwardx = isnan(e_forwardx) ? ZERO : e_forwardx;
      e_backwardx = isnan(e_backwardx) ? ZERO : e_backwardx;
      e_forwardy = isnan(e_forwardy) ? ZERO : e_forwardy;
      e_backwardy = isnan(e_backwardy) ? ZERO : e_backwardy;

      realtype extensible_x = params.gamma*(e_forwardx + e_backwardx);
      realtype extensible_y = params.gamma*(e_forwardy + e_backwardy);

      // Compute lower-substrate Pressure
      realtype p1, p2, p4, p5, p11;
      p1 = params.sigma / y;
      p2 = p1*p1;
      p4 = p2*p2;
      p5 = p4*p1;
      p11 = p5*p5*p1;

      realtype vdW_y = -(params.pressure*PI*(TWO*p11-FOUR*p5));

      // Cummulative fiber-dependent force
      total_fiber_force_x = -(bending_x + extensible_x);
      total_fiber_force_y = -(bending_y + extensible_y + vdW_y);

      NV_Ith_S(udot, xindex) += total_fiber_force_x;
      NV_Ith_S(udot, yindex) += total_fiber_force_y;
    }
  }

  return 0;
}

void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd_fiber,
                    std::vector<int>& nhbd_partner,
                    std::vector<bool>& mask,
                    parameter& param)
{
  nhbd_fiber.clear();
  nhbd_partner.clear();
  mask.clear();

  realtype r2 = param.radius * param.radius;
  int size = param.n * param.m;

  realtype x1, x2, y1, y2;
  for(int i = 0; i < size; ++i)
  {
    for(int j = i+1; j < size; ++j)
    {
      int fiber_i = i / param.n;
      int fiber_j = j / param.n;
      int pos_i = i % param.n;
      int pos_j = j % param.n;

      if(fiber_i != fiber_j
        or (pos_i != pos_j and pos_i != pos_j+1 and pos_i != pos_j-1))
      {
        x1 = NV_Ith_S(u, i);
        y1 = NV_Ith_S(u, i + size);
        x2 = NV_Ith_S(u, j);
        y2 = NV_Ith_S(u, j + size);

        realtype dist2 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
        if(dist2 <= r2)
        {
          nhbd_fiber.push_back(i);
          nhbd_partner.push_back(j);
          mask.push_back(true);
          mask.push_back(true);
        }
      }
    }
  }

  for(int i = 0; i < size; ++i)
  {
    for(int j = 0; j < param.sub_count; ++j)
    {
      x1 = NV_Ith_S(u, i);
      y1 = NV_Ith_S(u, i + size);
      x2 = NV_Ith_S(u, 2*size) + j*param.sub_h;
      y2 = NV_Ith_S(u, 2*size + 1);

      realtype dist2 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
      if(dist2 <= r2)
      {
        nhbd_fiber.push_back(i);
        nhbd_partner.push_back(j);
        mask.push_back(true);
        mask.push_back(false);
      }
    }
  }

  for(int i = 0; i < size; ++i)
  {
    for(int j = 0; j < param.osub_count; ++j)
    {
      x1 = NV_Ith_S(u, i);
      y1 = NV_Ith_S(u, i + size);
      x2 = param.osub + j*param.osub_h;
      y2 = ZERO;

      realtype dist2 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
      if(dist2 <= r2)
      {
        nhbd_fiber.push_back(i);
        nhbd_partner.push_back(j);
        mask.push_back(false);
        mask.push_back(true);
      }
    }
  }
}
