#ifndef SHEET1_H
#define SHEET1_H

#include "../algo/convergence.h"
#include "../base/types.h"
#include "../ivp/euler.h"
#include "test_runge_kutta.h"

template <typename Vector>
Vector RHS_P11(typename Vector::value_type t, const Vector &u)
{
  return u;
}

template <typename Vector>
Vector RHS_P13(typename Vector::value_type t, const Vector &u)
{
  return t * u;
};

void Problem_P12(FP_Type h)
{
  // IVP y'(t) = c*y(t) on the interval [0, 2]
  std_tWrapper f(RHS_P11<VectorD2>, 1);
  FP_Type t0 = 0.0;
  FP_Type t1 = 2.0;

  // initial value u(0) = 1
  dealii::Vector<FP_Type> u0(1);
  u0[0] = 1.0;

  // exact solution of the IVP at t = 2
  dealii::Vector<FP_Type> u(1);
  u[0] = std::exp(2);

  Blackbox Method(f, t0, u0);
  FP_Type EOC = eoc_1step(Method, t1, h);
  FP_Type OOC = ooc_1step(Method, t1, h, u);

  std::cout << "EOC: (Blackbox, h = " << h << ") "
            << EOC << std::endl
            << "OOC: (Blackbox, h = " << h << ") "
            << OOC << std::endl;

  ERK<DOPRI> Method2(f, t0, u0);
  FP_Type EOC2 = eoc_1step(Method2, t1, h);
  FP_Type OOC2 = ooc_1step(Method2, t1, h, u);

  std::cout << "EOC: (DOPRI, h = " << h << ") "
            << EOC2 << std::endl
            << "OOC: (DOPRI, h = " << h << ") "
            << OOC2 << std::endl;
}

void Problem_P13(FP_Type h)
{
  // IVP u'(t) = t * u(t) on the interval [0, 1]
  FP_Type t0 = 0.0;
  FP_Type t1 = 1.0;

  // Initial value u(0) = PI
  dealii::Vector<FP_Type> u0(1);
  u0[0] = M_PI;

  // Exact solution of the IVP at t = 1
  dealii::Vector<FP_Type> u(1);
  u[0] = M_PI * std::sqrt(std::exp(1));

  std_tWrapper f(RHS_P13<VectorD2>, 1);
  Euler Method(f, t0, u0);
  Method.iterate_up_to(t1, h);
  size_t steps = Method.n_steps();

  // Accuracy of computed approximation
  dealii::Vector<FP_Type> y = Method.approx();
  FP_Type diff = (y - u).l2_norm();

  std::cout << "Evaluations: (Euler, h = " << h << ") "
            << steps << std::endl
            << "Norm: |y - u| = "
            << diff << std::endl;
}

void Test_Sheet1()
{
  // Too numerically unstable for smaller values of h.
  Problem_P12(1e-1);
  std::cout << std::endl;

  Problem_P13(1e-1);
  Problem_P13(1e-2);
  Problem_P13(1e-3);
}

#endif // SHEET1_H
