#ifndef SHEET1_H
#define SHEET1_H
#include "../base/number_type.h"
#include "../ivp/blackbox.h"
#include "../ivp/eoc.h"

void Problem_P11(NumberType h)
{
  // IVP y'(t) = c*y(t) on the interval [0, 2]
  NumberType t0 = 0.0;
  NumberType t1 = 2.0;
  RHS_P11 rhs;

  // initial value u(0) = 1
  dealii::Vector<NumberType> u0(1);
  u0[0] = 1.0;

  // exact solution of the IVP at t = 2
  dealii::Vector<NumberType> u1(1);
  u1[0] = std::exp(2);

  evaluate_with_eoc<Blackbox>(rhs, t0, t1, h, u0, u1);
}

void Problem_P13(NumberType h)
{
  // IVP u'(t) = t * u(t) on the interval [0, 1]
  NumberType t0 = 0.0;
  NumberType t1 = 1.0;
  RHS_P13 rhs;

  // initial value u(0) = pi
  dealii::Vector<NumberType> u0(1);
  u0[0] = M_PI;

  // exact solution of the IVP at t = 1
  dealii::Vector<NumberType> u1(1);
  u1[0] = M_PI * std::sqrt(std::exp(1));

  evaluate_with_eoc<Explicit_Euler>(rhs, t0, t1, h, u0, u1);
}

void run_sheet1()
{
  std::cout << "Blackbox method (P11), h = 1e-2" << std::endl;
  Problem_P11(1e-2);
  std::cout << std::endl << "Blackbox method (P11), h = 1e-3" << std::endl;
  Problem_P11(1e-3);
  std::cout << std::endl << "Blackbox method (P11), h = 1e-4" << std::endl;
  Problem_P11(1e-4);

  std::cout << std::endl << "Euler method (P13), h = 1e-2" << std::endl;
  Problem_P13(1e-2);
  std::cout << std::endl << "Euler method (P13), h = 1e-3" << std::endl;
  Problem_P13(1e-3);
  std::cout << std::endl << "Euler method (P13), h = 1e-4" << std::endl;
  Problem_P13(1e-4);
}

#endif // SHEET1_H
