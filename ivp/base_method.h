#ifndef IVP_METHOD_H
#define IVP_METHOD_H
#include <cassert>
#include <vector>
#include <iostream>
#include <deal.II/base/data_out_base.h>
#include "../base/number_type.h"
#include "../base/functor.h"

/* This class solves an IVP of shape
 *     u'(t) = f(t, u(t));  u(t_0) = u_0
 *
 * where u_0.size() determines the problem dimension.
 *
 * Constructor:
 *    f       (functor) Right hand side of ODE in standard form.
 *    t0      (NumberType)  Initial time value. Default is 0.
 *    y0      (vector)  Initial value. Default is 1.
 *    nsteps  (int)     NumberType of integration steps. Default is 100.
 *    h       (NumberType)  Step length. Default is 1e-2.
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
  std::vector<NumberType> timepoints;

  // Use dealii vector for numerical operations
  std::vector<dealii::Vector<NumberType> > uapprox;

  // Constructor with initial value arguments
  IVP_Method(Functor &_f, NumberType _t0, dealii::Vector<NumberType> _u0, NumberType _h = 1.0e-2) :
    f(_f), t0(_t0), u0(_u0), h(_h)
  {}

  // Constructor with default initial values (t0 = 0, u0 = [1])
  IVP_Method(Functor &_f, NumberType _h = 1.0e-2) :
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

  virtual void iteration_step(dealii::Vector<NumberType> &u, NumberType &t, const NumberType &h)
  {
    throw std::invalid_argument("Please specify the step procedure in a child class");
  }

  // Execute the iteration over the given time interval [t0, t1], where t0
  // is set in the constructor, and t1 (t_limit) as argument to this function.
  void iterate_up_to(NumberType t_limit, bool adaptive_length = false)
  {
    // init iteration variables
    NumberType t = t0;
    NumberType steps_quotient = (t_limit - t0) / h;
    dealii::Vector<NumberType> u = u0;

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
    NumberType remainder = t_limit - t;
    if (remainder > 0)
      {
        nsteps++;
        iteration_step(u, t, remainder);

        timepoints.push_back(t);
        uapprox.push_back(u);
      }
  }

private:
  Functor &f;
  NumberType t0;
  dealii::Vector<NumberType> u0;
  NumberType h;
  size_t nsteps;
};

#endif // IVP_METHOD_H
