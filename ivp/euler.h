#ifndef EXPLICIT_EULER_H
#define EXPLICIT_EULER_H
#include <iostream>

#include "eos_method.h"
#include "../base/types.h"

class Explicit_Euler : public EOS_Method
{
  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using EOS_Method::EOS_Method;

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    return f(t, u);
  }

private:
  // EOS_Method::f;
  // EOS_Method::t0;
  // EOS_Method::u0;
  // EOS_Method::timepoints;
  // EOS_Method::uapprox;
  // EOS_Method::h;
  // EOS_Method::steps;
};

#endif // EXPLICIT_EULER_H
