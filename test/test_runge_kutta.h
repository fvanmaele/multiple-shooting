#ifndef TEST_RUNGE_KUTTA_H
#define TEST_RUNGE_KUTTA_H

#include "../lac/vector_operators.h"
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

    return 1./6*k1 + 2./6*k2 + 2./6*k3 + 1./6*k4;
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

class Blackbox : public EOS_Method
{
public:
  // Coefficients for intermediate computations and iteration steps
  // NOTE: With C++17, "static inline const" may be used
  const std::vector<FP_Type> c2 = { 1./5, 1./5 };
  const std::vector<FP_Type> c3 = { 3./10, 3./40, 9./40 };
  const std::vector<FP_Type> c4 = { 4./5, 44./45, 56./15, 32./9 };
  const std::vector<FP_Type> c5 = { 8./9, 19372./6561, 25360./2187,
                                    64448./6561, 212./729 };
  const std::vector<FP_Type> c6 = { 1., 9017./3168, 355./33,
                                    46732./5247, 49./176, 5103./18656 };

  const FP_Type s1 = 35./384;
  const FP_Type s2 = 0.;
  const FP_Type s3 = 500./1113;
  const FP_Type s4 = 125./192;
  const FP_Type s5 = 2187./6784;
  const FP_Type s6 = 11./84;

  // For base-class constructors, C++11 allows a class to specify that
  // base class constructors will be inherited.
  using EOS_Method::EOS_Method;

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    // intermediate computations for current step
    dealii::Vector<FP_Type> k1, k2, k3, k4, k5, k6;
    // functor return type: dealii::Vector
    k1 = f(t, u);
    k2 = f(t + c2[0]*h, u + h*(c2[1]*k1));
    k3 = f(t + c3[0]*h, u + h*(c3[1]*k1 + c3[2]*k2));
    k4 = f(t + c4[0]*h, u + h*(c4[1]*k1 - c4[2]*k2 + c4[3]*k3));
    k5 = f(t + c5[0]*h, u + h*(c5[1]*k1 - c5[2]*k2 + c5[3]*k3 - c5[4]*k4));
    k6 = f(t + c6[0]*h, u + h*(c6[1]*k1 - c6[2]*k2 + c6[3]*k3 + c6[4]*k4 - c6[5]*k5));

    // return increment
    return s1*k1 + s3*k3 + s4*k4 - s5*k5 + s6*k6;
  }

private:
  // EOS_Method::f;
  // EOS_Method::t0;
  // EOS_Method::u0;
  // EOS_Method::steps;
  // EOS_Method::timepoints;
  // EOS_Method::uapprox;
};

#endif // TEST_RUNGE_KUTTA_H
