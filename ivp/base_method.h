#ifndef IVP_METHOD_H
#define IVP_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>

#include <deal.II/base/data_out_base.h>

#include "../base/functor.h"
#include "../base/types.h"

/* This class solves an IVP of shape
 *     u'(t) = f(t, u(t));  u(t_0) = u_0
 *
 * where u_0.size() determines the problem dimension.
 *
 * Constructor:
 *    f       (functor) Right hand side of ODE in standard form.
 *    t0      (FP_Type)  Initial time value. Default is 0.
 *    y0      (vector)  Initial value. Default is 1.
 *    nsteps  (int)     FP_Type of integration steps. Default is 100.
 *    h       (FP_Type)  Step length. Default is 1e-2.
 *
 * The common wrapped functionality includes collecting of intermediary
 * computation results.
*/
class IVP_Method
{
public:
  // Implementations with fixed step length (template method pattern)
  friend class Blackbox;
  friend class Explicit_Euler;
  friend class ERK;
  friend class ERK_Test_O4;

  // Result vectors
  std::vector<FP_Type> timepoints;

  // Use dealii vector for numerical operations
  std::vector<dealii::Vector<FP_Type> > uapprox;

  // Constructor with initial value arguments
  IVP_Method(RHS &_f, FP_Type _t0, dealii::Vector<FP_Type> _u0, FP_Type _h = 1.0e-2) :
    f(_f), t0(_t0), u0(_u0), h(_h)
  {}

  // Constructor with default initial values (t0 = 0, u0 = [1])
  IVP_Method(RHS &_f, FP_Type _h = 1.0e-2) :
    f(_f), t0(0.0), h(_h)
  {
    u0.reinit(1);
    u0[0] = 1.0;
  }

  size_t steps() const
  {
    return nsteps;
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

  virtual void iteration_step(dealii::Vector<FP_Type> &u, FP_Type &t, const FP_Type &h)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  // Execute the iteration over the given time interval [t0, t1], where t0
  // is set in the constructor, and t1 (t_limit) as argument to this function.
  void iterate_up_to(FP_Type t_limit, bool adaptive_length = false)
  {
    // init iteration variables
    FP_Type t = t0;
    FP_Type steps_quotient = (t_limit - t0) / h;
    dealii::Vector<FP_Type> u = u0;

    // Use unsigned int to allow for higher amount of steps; for safe use,
    // add a check that the quotient I/h is not < 0.
    if (steps_quotient < 0)
        throw std::invalid_argument("Invalid step amount (are [t0, t1], h correct?)");
    nsteps = static_cast<size_t>(steps_quotient);

    // init output variables
    timepoints.assign(1, t);
    uapprox.assign(1, u);

    // continue iteration until either the amount of steps is reached,
    // or the iteration step would exceed the time interval.
    for (size_t i = 1; i <= nsteps; i++)
      {
        // write new values for u and t depending on iteration step
        iteration_step(u, t, h);

        // save values accordingly
        timepoints.push_back(t);
        uapprox.push_back(u);
      }

    // guarantee to exactly hit the right interval end
    FP_Type remainder = t_limit - t;
    if (remainder > 0)
      {
        nsteps++;
        iteration_step(u, t, remainder);

        timepoints.push_back(t);
        uapprox.push_back(u);
      }
  }

private:
  RHS &f;
  FP_Type t0;
  dealii::Vector<FP_Type> u0;

  FP_Type h;
  size_t nsteps;
};

#endif // IVP_METHOD_H
