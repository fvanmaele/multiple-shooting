#ifndef SHEET2_H
#define SHEET2_H

#include <cassert>

#include "../algo/convergence.h"
#include "../base/functor.h"
#include "../base/types.h"
#include "../ivp/euler.h"
#include "../ivp/runge_kutta.h"

void Problem_P21(FP_Type h, std::ostream &output)
{
  RHS_P21 rhs(5, 2, 2, 1);
  FP_Type t0 = 0.0;
  FP_Type t1 = 10.0;

  dealii::Vector<FP_Type> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the explicit Euler method
  Euler Method1(rhs, t0, u0);
  Method1.iterate_up_to(t1, h);
  Method1.print(output);
}

void Problem_P22(FP_Type h, std::ostream &output)
{
  RHS_P21 rhs(5, 2, 2, 1);
  FP_Type t0 = 0.0;
  FP_Type t1 = 10.0;

  dealii::Vector<FP_Type> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the classic Runge-Kutta method
  ERK_04 Tableau;
  dealii::LAPACKFullMatrix<FP_Type> A = MFA<16>(4, 4, Tableau.A);

  dealii::Vector<FP_Type> b(Tableau.b.begin(), Tableau.b.end());
  dealii::Vector<FP_Type> c(Tableau.c.begin(), Tableau.c.end());

  ERK Method2(rhs, t0, u0, A, b, c);
  Method2.iterate_up_to(t1, h);
  Method2.print(output);

  // The classical Range-Kutta method has order of convergence 4.
  FP_Type EOC = eoc_1step(Method2, t1, 1e-2, 2);
  std::cout << "EOC: (Range-Kutta, h = " << 1e-2 << ") "
            << EOC << std::endl;
}

void Test_Sheet2()
{
  std::cout << std::endl
            << "Lotka-Volterra equations (P21), h = 0.1" << std::endl;

  std::ofstream output_file;
  output_file.open("volterra_euler.dat");
  assert(output_file.is_open());

  // Generate output data and run GNUPLOT (via shell)
  Problem_P21(1e-1, output_file);
  output_file.close();
  std::system("gnuplot -p -e \"plot 'volterra_euler.dat' using 1:2 with lines, "
              "'volterra_euler.dat' using 1:3 with lines\"");

  // Generate output data for Runge-Kutta method
  output_file.open("volterra_runge.dat");
  assert(output_file.is_open());

  Problem_P22(1e-1, output_file);
  output_file.close();
  std::system("gnuplot -p -e \"plot 'volterra_runge.dat' using 1:2 with lines, "
              "'volterra_runge.dat' using 1:3 with lines\"");
}

#endif // SHEET2_H
