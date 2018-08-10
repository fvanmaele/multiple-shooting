#ifndef EOS_METHOD_H
#define EOS_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>

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
  OneStepMethod(TimeFunctor &_f, FP_Type _t0, VectorD2 _u0, Curve *_u = nullptr)
    :
      f(_f), u(_u), t0(_t0), u0(_u0), steps(0), Yn(_u0.size()),
      timepoints(1, _t0), uapprox(1, _u0)
  {}

  // Access functions
  VectorD2 approx() const
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

  MatrixD2 fund_matrix() const
  {
    return Yn;
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

  bool sol_is_nan(const VectorD2 &y)
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
    steps = 0;
  }

  void save_step(const FP_Type &t, const VectorD2 &u)
  {
    timepoints.push_back(t);
    uapprox.push_back(u);
  }

  virtual VectorD2
  increment_function(FP_Type, const VectorD2&, FP_Type)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  virtual std::pair<VectorD2, MatrixD2>
  increment_variational(FP_Type, const VectorD2&, FP_Type, const MatrixD2&, TimeDivFunctor*)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  virtual ~OneStepMethod() = default;

  // Execute the iteration over the given time interval [t0, t_limit].
  void iterate_up_to(FP_Type t_lim, FP_Type h,
                     bool fundamental_matrix = false,
                     FP_Type C = 2)
  {
    reset(); // init output variables

    FP_Type  t = t0;
    VectorD2 y = u0;
    size_t n = u0.size();

    // Check prerequisites for variational equation
    TimeDivFunctor* f_diff = dynamic_cast<TimeDivFunctor*>(&f);

    if (fundamental_matrix && f_diff == nullptr)
      throw std::invalid_argument("right-hand side is not differentiable");

    // Initial value for Yn(t; t0, u0)
    Yn = dealii::IdentityMatrix(n);

    while (t_lim - t > 0)
      { // Avoid rounding errors (Rem. 2.4.3)
        if (t + 1.1*h >= t_lim)
          h = t_lim - t;

        if (fundamental_matrix)
          {
            auto U = increment_variational(t, y, h, Yn, f_diff);
            y += h * U.first;
            Yn.add(1, h * U.second);
          }
        else
          { // y_k = y_{k-1} + h*F(t_{k-1}, y_{k-1})
            y += h * increment_function(t, y, h);
          }
        // t_k = t_{k-1} + h
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

    if (timepoints.back() != t_lim)
      {
        std::string err = "time step outside interval end ("
            + std::to_string(timepoints.back()) + "; ["
            + std::to_string(t0) + ", "
            + std::to_string(t_lim) + "])";

        throw std::out_of_range(err.c_str());
      }
  }

private:
  TimeFunctor &f;  // rhs of ODE in standard form
  Curve *u;        // exact solution
  FP_Type t0;      // initial time value
  VectorD2 u0;     // initial value
  size_t steps;    // amount of integration steps
  MatrixD2 Yn;     // fundamental matrix

  // Result vectors for the IVP
  std::vector<FP_Type> timepoints;
  std::vector<VectorD2> uapprox;
};

#endif // EOS_METHOD_H
