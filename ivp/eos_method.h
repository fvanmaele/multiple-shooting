#ifndef EOS_METHOD_H
#define EOS_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>

#include <deal.II/base/data_out_base.h>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"

/* This class solves an IVP of shape
 *     u'(t) = f(t, u(t));  u(t_0) = u_0
 *
 * using an explicit one-step method.
 *
 * The common wrapped functionality includes collecting of intermediary
 * computation results.
*/
class OneStepMethod
{
public:
  // Implementations with fixed step length (template method pattern)
  friend class Blackbox;
  friend class Euler;
  friend class ERK_Test_04;

  template <typename BTab>
  friend class ERK;

  // The exact solution may be specified as an optional argument
  // for comparison purposes.
  OneStepMethod(TimeFunctor &_f, FP_Type _t0, dealii::Vector<FP_Type> _u0,
             Curve *_u = nullptr)
    :
      f(_f), u(_u), t0(_t0), u0(_u0),
      timepoints(1, _t0), uapprox(1, _u0)
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

  bool sol_is_nan(const dealii::Vector<FP_Type> &y)
  {
    for (size_t i = 0; i < y.size(); i++)
      {
        if (std::isnan(y[i]))
          {
            return true;
          }
      }
    return false;
  }

  void reset()
  {
    timepoints.assign(1, t0);
    uapprox.assign(1, u0);
  }

  void save_step(const FP_Type &t, const dealii::Vector<FP_Type> &u)
  {
    timepoints.push_back(t);
    uapprox.push_back(u);
  }

  virtual dealii::Vector<FP_Type>
  increment_function(FP_Type t, const dealii::Vector<FP_Type> &y, FP_Type h)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  virtual ~OneStepMethod() = default;

  // Execute the iteration over the given time interval [t0, t_limit].
  void iterate_up_to(FP_Type t_lim, FP_Type h, FP_Type C = 2)
  {
    // init input variables
    FP_Type t = t0;
    FP_Type h_arg = h; // copy for step check
    dealii::Vector<FP_Type> y = u0;

    // init output variables
    reset();
    steps = 0;

    while (t_lim - t > 0)
      { // Avoid rounding errors (Rem. 2.4.3)
        if (t + 1.1*h >= t_lim)
          h = t_lim - t;

        // y_k = y_{k-1} + h*F(t_{k-1}, y_{k-1})
        // t_k = t_{k-1} + h
        y += h * increment_function(t, y, h);
        t += h;

        // Add u_k, t_k to result vectors
        save_step(t, y);
        steps++;

        if (sol_is_nan(y))
          throw std::overflow_error("global error too large (NaN)");

        if (u != nullptr)
          if (y.l2_norm() >= C*(*u)(t).l2_norm())
            throw std::domain_error("global error too large");
      }

    if (static_cast<FP_Type>(steps) != (t_lim - t0) / h_arg)
      std::cerr << "warning: mismatch in step amount ("
                << std::setprecision(15) << std::setw(15)
                << (t_lim - t0) / h_arg << " to " << steps << ")"
                << std::endl;

    if (timepoints.back() != t_lim)
      std::cerr << "warning: time step outside interval end ("
                << timepoints.back() << "; [" << t0 << ", " << t_lim << "])"
                << std::endl;
  }

private:
  TimeFunctor &f;             // rhs of ODE in standard form
  Curve *u;                   // exact solution
  FP_Type t0;                 // initial time value
  dealii::Vector<FP_Type> u0; // initial value
  size_t steps;               // amount of integration steps

  // Result vectors for the IVP
  std::vector<FP_Type> timepoints;
  std::vector<dealii::Vector<FP_Type> > uapprox;
};

#endif // EOS_METHOD_H
