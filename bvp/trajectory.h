#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>

#include <deal.II/base/function_cspline.h>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../ivp/runge_kutta.h"

/*!
 * \brief Find suitable initial values for subintervals, using an approximate
 * solution to the BVP (a \e starting \e trajectory).
 */
template <typename M = DOPRI54>
std::vector<FP_Type>
trajectory(FP_Type a, FP_Type b, TimeFunctor &rhs, Curve *eta,
           FP_Type C = 2, bool ssc = true,
           FP_Type h0 = 1e-1, FP_Type TOL = 1e-4)
{
  assert(a < b);
  std::vector<FP_Type> T;

  FP_Type t_i = a;
  T.push_back(t_i);

  while (t_i < b)
    {
      VectorD2 eta_i = (*eta)(t_i);
      ERK<M> AdM(rhs, t_i, eta_i, false, eta);

      try
      {
        if (ssc)
          AdM.iterate_with_ssc(b, h0, TOL, C);
        else
          AdM.iterate_up_to(b, h0, C);

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

/*!
 * \brief Create a vector of equally spaced points.
 *
 * Based on MATLAB linspace.
 */
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

/*!
 * \brief Interpolate a set of intervals to a given amount \f$N\f$
 * using Spline interpolation.
 *
 * Requires GSL support.
 */
std::vector<FP_Type>
interpolate_points(const std::vector<FP_Type> &t, size_t N)
{
  const size_t k = t.size();
  std::vector<FP_Type> points(k);
  std::vector<FP_Type> values(k);

  for (size_t i = 0; i < k; i++)
    {
      points[i] = i/(FP_Type)(k-1);
      values[i] = t[i];
    }

  dealii::Functions::CSpline<1> S(points, values);
  std::vector<FP_Type> result(N);
  std::vector<FP_Type> u = linspace(0, 1, N);

  for (size_t i = 0; i < u.size(); i++)
    {
      dealii::Point<1, FP_Type> x(u.at(i));
      result.at(i) = S.value(x);
    }
  return result;
}

#endif // TRAJECTORY_H
