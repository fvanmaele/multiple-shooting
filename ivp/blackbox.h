#ifndef BLACKBOX_METHOD_H
#define BLACKBOX_METHOD_H
#include <vector>
#include "../base/functor.h"
#include "../base/number_type.h"
#include "../ivp/base_method.h"
#include "../lac/vector_operators.h"

class Blackbox : public IVP_Method
{
public:
  // Coefficients for intermediate computations and iteration steps
  // NOTE: With C++17, "static inline const" may be used
  const std::vector<NumberType> c2 = { 1./5, 1./5 };
  const std::vector<NumberType> c3 = { 3./10, 3./40, 9./40 };
  const std::vector<NumberType> c4 = { 4./5, 44./45, 56./15, 32./9 };
  const std::vector<NumberType> c5 = { 8./9, 19372./6561, 25360./2187,
                                   64448./6561, 212./729 };
  const std::vector<NumberType> c6 = { 1., 9017./3168, 355./33,
                                          46732./5247, 49./176, 5103./18656 };

  const NumberType s1 = 35./384;
  const NumberType s2 = 0.;
  const NumberType s3 = 500./1113;
  const NumberType s4 = 125./192;
  const NumberType s5 = 2187./6784;
  const NumberType s6 = 11./84;

  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using IVP_Method::IVP_Method;

  virtual void iteration_step(dealii::Vector<NumberType> &u, NumberType &t, const NumberType &h) override
  {
    // intermediate computations for current step
    dealii::Vector<NumberType> k1, k2, k3, k4, k5, k6;
    k1 = f(t, u); // functor return type: dealii::Vector
    k2 = f(t + c2[0]*h, u + h*(c2[1]*k1));
    k3 = f(t + c3[0]*h, u + h*(c3[1]*k1 + c3[2]*k2));
    k4 = f(t + c4[0]*h, u + h*(c4[1]*k1 - c4[2]*k2 + c4[3]*k3));
    k5 = f(t + c5[0]*h, u + h*(c5[1]*k1 - c5[2]*k2 + c5[3]*k3 - c5[4]*k4));
    k6 = f(t + c6[0]*h, u + h*(c6[1]*k1 - c6[2]*k2 + c6[3]*k3 + c6[4]*k4 - c6[4]*k5));

    // perform step
    u = u + h*(s1*k1 + s3*k3 + s4*k4 - s5*k5 + s6*k6);
    t = t + h;
  }

private:
  // IVP_Method::f;
  // IVP_Method::t0;
  // IVP_Method::u0;
  // IVP_Method::h;
  // IVP_Method::nsteps;
};

#endif // BLACKBOX_METHOD_H
