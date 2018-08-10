#ifndef TEST_VAREQ_H
#define TEST_VAREQ_H

#include <cmath>

#include "../base/forward_ad.h"
#include "../ivp/runge_kutta.h"

#define VAREQ_METHOD DOPRI54

namespace Test
{
  VectorAD RHS_Var1(NumberAD, const VectorAD &u)
  {
    VectorAD y(1);
    y[0] = -std::exp(u[0]);
    return y;
  }

  VectorAD RHS_Var2(NumberAD, const VectorAD &u)
  {
    VectorAD y(1);
    y[0] = -std::pow(u[0], 2);
    return y;
  }

  VectorAD RHS_Var3(NumberAD t, const VectorAD &u)
  {
    VectorAD y(1);
    y[0] = u[0] + std::sqrt(t*t + u[0]*u[0]) / t;
    return y;
  }

  VectorAD RHS_Var4(NumberAD, const VectorAD &u)
  {
    VectorAD y(1);
    y[0] = -0.5 * std::pow(u[0], 3);
    return y;
  }

  void Var1()
  {
    FAD_tWrapper f_ad(RHS_Var1, 1);
    FP_Type TOL = 1e-8;
    FP_Type t0 = 1;
    FP_Type t1 = 20;

    VectorD2 u0(1);
    u0[0] = 0;

    // Exact solution - IVP
    VectorD2 u(1);
    u[0] = -std::log(t1);

    // Exact solution - Variational equation
    MatrixD2 phi(1, 1);
    phi(0, 0) = 1./t1;

    ERK<VAREQ_METHOD> Method(f_ad, t0, u0);
    Method.iterate_with_ssc(t1, 1e-1, TOL, true);

    VectorD2 y = Method.approx();
    MatrixD2 phi_y = Method.fund_matrix();

    std::cout << "Exact solution (IVP): "
              << std::endl << u
              << "Approximate: "
              << std::endl << y
              << "Exact solution (VarEq): "
              << std::endl << phi
              << "Approximate: "
              << std::endl << phi_y;
  }

  void Var2()
  {
    FAD_tWrapper f_ad(RHS_Var2, 1);
    FP_Type TOL = 1e-8;
    FP_Type t0 = 1;
    FP_Type t1 = 20;

    VectorD2 u0(1);
    u0[0] = 1;

    // Exact solution - IVP
    VectorD2 u(1);
    u[0] = 1./(t1+1);

    // Exact solution - Variational equation
    MatrixD2 phi(1, 1);
    phi(0, 0) = 1./std::pow(t1+1, 2);

    ERK<VAREQ_METHOD> Method(f_ad, t0, u0);
    Method.iterate_with_ssc(t1, 1e-1, TOL, true);

    VectorD2 y = Method.approx();
    MatrixD2 phi_y = Method.fund_matrix();

    std::cout << "Exact solution (IVP): "
              << std::endl << u
              << "Approximate: "
              << std::endl << y
              << "Exact solution (VarEq): "
              << std::endl << phi
              << "Approximate: "
              << std::endl << phi_y;
  }

  void Var3()
  {
    FAD_tWrapper f_ad(RHS_Var3, 1);
    FP_Type TOL = 1e-8;
    FP_Type t0 = 1;
    FP_Type t1 = 20;

    VectorD2 u0(1);
    u0[0] = 0;

    // Exact solution - IVP
    VectorD2 u(1);
    u[0] = 0.5 * (t1*t1 - 1);

    // Exact solution - Variational equation
    MatrixD2 phi(1, 1);
    phi(0, 0) = t1*t1 + 1;

    ERK<VAREQ_METHOD> Method(f_ad, t0, u0);
    Method.iterate_with_ssc(t1, 1e-1, TOL, true);

    VectorD2 y = Method.approx();
    MatrixD2 phi_y = Method.fund_matrix();

    std::cout << "Exact solution (IVP): "
              << std::endl << u
              << "Approximate: "
              << std::endl << y
              << "Exact solution (VarEq): "
              << std::endl << phi
              << "Approximate: "
              << std::endl << phi_y;
  }

  void Var4()
  {
    FAD_tWrapper f_ad(RHS_Var4, 1);
    FP_Type TOL = 1e-8;
    FP_Type t0 = 1;
    FP_Type t1 = 20;

    VectorD2 u0(1);
    u0[0] = 1;

    // Exact solution - IVP
    VectorD2 u(1);
    u[0] = 1./std::sqrt(t1+1);

    // Exact solution - Variational equation
    MatrixD2 phi(1, 1);
    phi(0, 0) = 1./std::pow(t1+1, 1.5);

    ERK<VAREQ_METHOD> Method(f_ad, t0, u0);
    Method.iterate_with_ssc(t1, 1e-1, TOL, true);

    VectorD2 y = Method.approx();
    MatrixD2 phi_y = Method.fund_matrix();

    std::cout << "Exact solution (IVP): "
              << std::endl << u
              << "Approximate: "
              << std::endl << y
              << "Exact solution (VarEq): "
              << std::endl << phi
              << "Approximate: "
              << std::endl << phi_y;
  }
}

#endif // TEST_VAREQ_H
