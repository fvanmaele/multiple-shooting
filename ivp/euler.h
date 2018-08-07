#ifndef EULER_H
#define EULER_H
#include <iostream>

#include "eos_method.h"

class Euler : public OneStepMethod
{
  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using OneStepMethod::OneStepMethod;

  virtual dealii::Vector<FP_Type>
  increment_function(FP_Type t, const dealii::Vector<FP_Type> &y, FP_Type) override
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
