#ifndef SHEET1_H
#define SHEET1_H

#include "../algo/convergence.h"
#include "../base/types.h"
#include "../ivp/euler.h"
#include "test_runge_kutta.h"

namespace Test
{
  VectorD2 RHS_P11(FP_Type, const VectorD2 &u)
  {
    return u;
  }

  VectorD2 RHS_P13(FP_Type t, const VectorD2 &u)
  {
    return t * u;
  }

  void Problem_P12(FP_Type h)
  {
    // IVP y'(t) = c*y(t) on the interval [0, 2]
    std_tWrapper f(RHS_P11, 1);
    FP_Type t0 = 0.0;
    FP_Type t1 = 2.0;

    // initial value u(0) = 1
    VectorD2 u0(1);
    u0[0] = 1.0;

    // exact solution of the IVP at t = 2
    VectorD2 u(1);
    u[0] = std::exp(2);

    Blackbox M1(f, t0, u0);
    ConvergenceTester T1(M1, 3, 1e-1);

    std::cout << "EOC: (Blackbox, h = " << h << ") "
              << T1.eoc(t1, h) << std::endl
              << "OOC: (Blackbox, h = " << h << ") "
              << T1.ooc(t1, h, u) << std::endl;

    ERK<DOPRI54> M2(f, t0, u0);
    ConvergenceTester T2(M2, 3, 1e-1);
    
    std::cout << "EOC: (DOPRI54, h = " << h << ") "
              << T2.eoc(t1, h) << std::endl
              << "OOC: (DOPRI54, h = " << h << ") "
              << T2.ooc(t1, h, u) << std::endl;

    ERK<KARP> M3(f, t0, u0);
    ConvergenceTester T3(M3, 3, 1e-1);

    std::cout << "EOC: (KARP, h = " << h << ") "
              << T3.eoc(t1, h) << std::endl
              << "OOC: (KARP, h = " << h << ") "
              << T3.ooc(t1, h, u) << std::endl;

    ERK<ERK_04> M4(f, t0, u0);
    ConvergenceTester T4(M4, 3, 1e-1);

    std::cout << "EOC: (RK04, h = " << h << ") "
              << T4.eoc(t1, h) << std::endl
              << "OOC: (RK04, h = " << h << ") "
              << T4.ooc(t1, h, u) << std::endl;
  }

  void Problem_P13(FP_Type h)
  {
    // IVP u'(t) = t * u(t) on the interval [0, 1]
    FP_Type t0 = 0.0;
    FP_Type t1 = 1.0;

    // Initial value u(0) = PI
    VectorD2 u0(1);
    u0[0] = M_PI;

    // Exact solution of the IVP at t = 1
    VectorD2 u(1);
    u[0] = M_PI * std::sqrt(std::exp(1));

    std_tWrapper f(RHS_P13, 1);
    Euler Method(f, t0, u0);
    Method.iterate_up_to(t1, h);
    size_t steps = Method.n_steps();

    // Accuracy of computed approximation
    VectorD2 y = Method.approx();
    FP_Type diff = (y - u).l2_norm();

    std::cout << "Evaluations: (Euler, h = " << h << ") "
              << steps << std::endl
              << "Norm: |y - u| = "
              << diff << std::endl;
  }

  void Sheet1()
  {
    // Too numerically unstable for smaller values of h.
    Problem_P12(1e-1);
    Problem_P13(1e-1);
    Problem_P13(1e-2);
    Problem_P13(1e-3);
  }
}

#endif // SHEET1_H
