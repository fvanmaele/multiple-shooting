#ifndef EULER_H
#define EULER_H
#include <iostream>

#include "eos_method.h"
#include "../base/types.h"

class Euler : public EOS_Method
{
  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using EOS_Method::EOS_Method;

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    return f.value(t, u);
  }

private:
  // EOS_Method::f;
  // EOS_Method::t0;
  // EOS_Method::u0;
  // EOS_Method::steps;
  // EOS_Method::timepoints;
  // EOS_Method::uapprox;
  // EOS_Method::Y;
};

#endif // EULER_H
