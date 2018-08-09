#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../ivp/runge_kutta.h"

// Find suitable initial values for subintervals using an approximate
// solution to the BVP, i.e. a "starting trajectory".
std::vector<FP_Type>
trajectory(FP_Type a, FP_Type b, TimeFunctor &f, Curve *eta,
           FP_Type TOL = 1e-3, FP_Type C = 2, bool ssc = false)
{
  assert(a < b);
  std::vector<FP_Type> T;

  FP_Type t_i = a;
  T.push_back(t_i);

  while (t_i < b)
    {
      // Initial value of eta' = f(t, eta_i)
      VectorD2 eta_i = (*eta)(t_i);
      ERK<DOPRI87> AdM(f, t_i, eta_i, eta);

      try
      {
        if (ssc)
          AdM.iterate_with_ssc(b, 1e-3, TOL, false, C);
        else
          AdM.iterate_up_to(b, TOL, C);

        t_i = AdM.endpoint();
        t_i < b ? T.push_back(t_i) : T.push_back(b);
      }
      catch (std::domain_error &e)
      {
        t_i = AdM.endpoint();
        t_i < b ? T.push_back(t_i) : T.push_back(b);
      }
    }

  if (T.back() != b)
    T.push_back(b);

  return T;
}

// Create vector of equidistant points
std::vector<FP_Type>
linspace(FP_Type a, FP_Type b, int n = 100)
{
  assert(a < b);
  assert(n > 0);

  std::vector<FP_Type> v(n);
  FP_Type h = (b - a) / (FP_Type)(n-1);

  for (size_t i = 0; i < v.size(); i++)
    {
      v[i] = a + i*h;
    }
  return v;
}

#endif // TRAJECTORY_H
