#ifndef EOS_METHOD_H
#define EOS_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>

#include <deal.II/base/data_out_base.h>

#include "../base/types.h"
#include "../base/functor.h"
#include "../base/forward_ad.h"

/* This class solves an IVP of shape
 *     u'(t) = f(t, u(t));  u(t_0) = u_0
 *
 * using an explicit one-step method.
 *
 * Constructor:
 *    f       (functor) Right hand side of ODE in standard form.
 *    t0      (FP_Type)  Initial time value. Default is 0.
 *    y0      (vector)  Initial value. Default is 1.
 *    steps  (int)     FP_Type of integration steps. Default is 100.
 *    h       (FP_Type)  Step length. Default is 1e-2.
 *
 * The common wrapped functionality includes collecting of intermediary
 * computation results.
*/
class EOS_Method
{
public:
  // Implementations with fixed step length (template method pattern)
  friend class Blackbox;
  friend class Euler;
  friend class ERK_Test_O4;

  template <typename BTab>
  friend class ERK;

  // The intial values could be stored solely in the timepoints
  // resp. uapprox vector, but this is left out for simplicity.
  EOS_Method(TimeFunctor &_f, FP_Type _t0, dealii::Vector<FP_Type> _u0)
    :
      f(_f), t0(_t0), u0(_u0),
      timepoints(1, _t0), uapprox(1, _u0),
      Y(dealii::IdentityMatrix(_u0.size())) // Y(t_0) = I
  {}

  // Access functions
  dealii::Vector<FP_Type> approx() const
  {
    return uapprox.back();
  }

  FP_Type endpoint() const
  {
    return timepoints.back();
  }

  dealii::FullMatrix<FP_Type> fund_matrix() const
  {
    return Y;
  }

  size_t n_steps() const
  {
    return steps;
  }

  // Print approximation at each step in a tabular format
  void print(std::ostream &out = std::cout) const
  {
    assert(timepoints.size() == uapprox.size());

    for (size_t i = 0; i < timepoints.size(); i++)
      {
        out << timepoints[i] << "\t" << uapprox[i];
      }
  }

  void reset()
  {
    timepoints.assign(1, t0);
    uapprox.assign(1, u0);
    Y = 0;
  }

  void save_step(const FP_Type &t,
                 const dealii::Vector<FP_Type> &u)
  {
    timepoints.emplace_back(t);
    uapprox.emplace_back(u);
  }

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h, TimeFunctor &func)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  virtual ~EOS_Method() = default;

  // Execute the iteration over the given time interval [t0, t_limit].
  void iterate_up_to(FP_Type t_lim, FP_Type h,
                     bool fundamental_matrix = false)
  {
    // init input variables
    FP_Type t = t0;
    dealii::Vector<FP_Type> u = u0;

    // init output variables
    reset();
    steps = 0;

    while (t_lim - t > 0)
      {
        // Guarantee to hit right interval end
        if (t + 1.1*h >= t_lim)
          h = t_lim - t;

        // Take last values (t, u) to compute step in variational equation,
        // including the corresponding step size.
        if (fundamental_matrix)
          // Note: the new step is written in-place, unlike the IVP solution.
          Y += h * fund_matrix_increment(t, u, h, Y);

        // u_k = u_{k-1} + h*F(t_{k-1}, u_{k-1})
        // t_k = t_{k-1} + h
        u += h * increment_function(t, u, h);
        t += h;

        // Add u_k, t_k to result vectors
        save_step(t, u);
        steps++;
      }
    assert(timepoints.back() == t_lim);
    assert(steps = (size_t)(t_lim - t0) / h);
  }

  // Compute increment function of variational equation from IVP step.
  dealii::FullMatrix<FP_Type>
  fund_matrix_increment(const FP_Type t, const dealii::Vector<FP_Type> &u,
                        const FP_Type h, const dealii::FullMatrix<FP_Type> &Y)
  {
    dealii::Vector<FP_Type> phi(Y.m());
    dealii::FullMatrix<FP_Type> Y_inc = Y;

    for (size_t j = 0; j < Y.n(); j++)
      {
        // Compute values column by column
        for (size_t i = 0; i < Y.m(); i++)
          phi(i) = Y(i, j);

        // Construct RHS of variational equation
        FundMatrixFunctor F(f, phi);
        phi = increment_function(t, u, h, F);

        // Write back result
        for (size_t i = 0; i < Y.m(); i++)
          Y_inc(i, j) = phi(i);
      }
    return Y_inc;
  }

private:
  TimeFunctor &f;
  FP_Type t0;
  dealii::Vector<FP_Type> u0;
  size_t steps;

  // Result vectors for the IVP
  std::vector<FP_Type> timepoints;
  std::vector<dealii::Vector<FP_Type> > uapprox;

  // Variational equation
  dealii::FullMatrix<FP_Type> Y;
};

#endif // EOS_METHOD_H