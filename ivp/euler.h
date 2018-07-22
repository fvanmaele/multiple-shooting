#ifndef EXPLICIT_EULER_H
#define EXPLICIT_EULER_H
#include <iostream>

#include "base_method.h"
#include "../base/types.h"

class Explicit_Euler : public IVP_Method
{
  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using IVP_Method::IVP_Method;

  virtual void iteration_step(dealii::Vector<FP_Type> &u, FP_Type &t, const FP_Type &h) override
  {
    u += h * f(t,u);
    t += h;
  }

private:
  // IVP_Method::f;
  // IVP_Method::t0;
  // IVP_Method::u0;
  // IVP_Method::h;
  // IVP_Method::nsteps;
};

#endif // EXPLICIT_EULER_H
