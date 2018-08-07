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

  // XXX: rounding errors: check boundary up to TOL
  while (std::round((b - t_i)/TOL)*TOL > TOL)
    {
      // Initial value of eta' = f(t, eta_i)
      dealii::Vector<FP_Type> eta_i = (*eta)(t_i);
      ERK<RK65> AdM(f, t_i, eta_i, eta);

      try
      {
        if (ssc)
          AdM.iterate_with_ssc(b, 1e-1, TOL, false, C);
        else
          AdM.iterate_up_to(b, TOL, C);
      }
      catch (std::domain_error &e)
      {
        t_i = AdM.endpoint();
        T.push_back(t_i);
      }
    }

  if (T.back() != b)
    T.push_back(b);

  return T;
}

#endif // TRAJECTORY_H
