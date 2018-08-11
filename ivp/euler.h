#ifndef EULER_H
#define EULER_H
#include <iostream>

#include "eos_method.h"

/*! \class Euler
 * \brief Classic Euler method of order 1.
 */
class Euler : public OneStepMethod
{
  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using OneStepMethod::OneStepMethod;

  virtual VectorD2
  increment_function(FP_Type t, const VectorD2 &y, FP_Type) override
  {
    return f(t, y);
  }

private:
  // OneStepMethod::f;
  // OneStepMethod::t0;
  // OneStepMethod::u0;
  // OneStepMethod::steps;
  // OneStepMethod::timepoints;
  // OneStepMethod::uapprox;
};

#endif // EULER_H
