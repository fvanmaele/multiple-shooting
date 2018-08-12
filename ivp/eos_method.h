#ifndef EOS_METHOD_H
#define EOS_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"

/*!
 * \brief Solve an IVP of shape \f$u'(t) = f(t, u(t)), u(t_0) = u_0\f$
 * using a one-step method.
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

  /*!
   * \brief OneStepMethod
   * \param _f Right-hand side of the ODE.
   * \param _t0 Initial time value.
   * \param _u0 Initial value, \f$u(t_0) = u_0\f$.
   * \param _var_eq If true, solve the variational equation.
   * \param _u Exact solution of the IVP.
   *
   * Solving the variational equation is done \e simultaneously with solving the
   * IVP. In particular, the same step width is used for both problems. See the
   * \c ERK documentation for further discussion.
   *
   * A provided \c exact solution \f$u(t)\f$ may be used to verify the local
   * error at each time step. The involved constant are specified to the
   * iteration methods.
   */
  OneStepMethod(TimeFunctor &_f, FP_Type _t0, VectorD2 _u0,
                bool _var_eq = false, Curve *_u = nullptr)
    :
      f(_f), u(_u), t0(_t0), u0(_u0), steps(0), f_div(nullptr), var_eq(_var_eq),
      timepoints(1, _t0), uapprox(1, _u0)
  {
    if (var_eq)
      { // Check prerequisites for variational equation
        f_div = dynamic_cast<TimeDivFunctor*>(&f);

        if (f_div == nullptr)
          throw std::invalid_argument("right-hand side is not differentiable");

        // Initial value for Yn(t; t0, u0)
        Yn = dealii::IdentityMatrix(u0.size());
      }
  }

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

  /*!
   * \brief Print the approximate solution at each time step in a tabular format.
   * \param File stream to print to, for example stdout or an output file.
   */
  void print(std::ostream &out = std::cout) const
  {
    assert(timepoints.size() == uapprox.size());

    for (size_t i = 0; i < timepoints.size(); i++)
      out << timepoints[i] << "\t" << uapprox[i];
  }

  /*!
   * \brief Check if an a \dealii::Vector element is NaN (\c std::isnan)
   */
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
  increment_variational(FP_Type, const VectorD2&, FP_Type, const MatrixD2&)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  virtual ~OneStepMethod() = default;

  // Execute the iteration over the given time interval [t0, t_limit].
  void iterate_up_to(FP_Type t_lim, FP_Type h, FP_Type C = 2)
  {
    reset(); // init output variables

    FP_Type  t = t0;
    VectorD2 y = u0;

    while (t_lim - t > 0)
      { // Avoid rounding errors (Rem. 2.4.3)
        if (t + 1.1*h >= t_lim)
          h = t_lim - t;

        if (var_eq)
          {
            auto U = increment_variational(t, y, h, Yn);

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
          throw std::overflow_error("local error too large (NaN)");

        if (u != nullptr && y.l2_norm() >= C*(*u)(t).l2_norm())
          throw std::range_error("local error too large");
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
  // Initial value problem
  TimeFunctor &f;  // rhs of ODE in standard form
  Curve *u;        // exact solution
  FP_Type t0;      // initial time value
  VectorD2 u0;     // initial value
  size_t steps;    // amount of integration steps

  // Variational equation
  TimeDivFunctor* f_div;
  MatrixD2 Yn;
  bool var_eq;

  // Result vectors for the IVP
  std::vector<FP_Type> timepoints;
  std::vector<VectorD2> uapprox;
};

#endif // EOS_METHOD_H
