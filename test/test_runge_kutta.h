#ifndef TEST_RUNGE_KUTTA_H
#define TEST_RUNGE_KUTTA_H

#include "../ivp/runge_kutta.h"

class ERK_Test_O4 : public EOS_Method
{
public:
  using EOS_Method::EOS_Method;

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    dealii::Vector<FP_Type> k1 = f(t, u);
    dealii::Vector<FP_Type> k2 = f(t + 0.5*h, u + 0.5*h*k1);
    dealii::Vector<FP_Type> k3 = f(t + 0.5*h, u + 0.5*h*k2);
    dealii::Vector<FP_Type> k4 = f(t + h, u + h*k3);

    return h * (1./6*k1 + 2./6*k2 + 2./6*k3 + 1./6*k4);
  }

private:
  // EOS_Method::f;
  // EOS_Method::t0;
  // EOS_Method::u0;
  // EOS_Method::timepoints;
  // EOS_Method::uapprox;
  // EOS_Method::h;
  // EOS_Method::nsteps;
};

#endif // TEST_RUNGE_KUTTA_H
