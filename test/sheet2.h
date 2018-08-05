#ifndef SHEET2_H
#define SHEET2_H

#include <cassert>

#include "../algo/convergence.h"
#include "../base/types.h"
#include "../ivp/euler.h"
#include "../ivp/runge_kutta.h"

VectorD2 LotkaVolterra(FP_Type, const VectorD2 &u)
{
  VectorD2 result(2);
  result[0] = u[0] * (5 - 2 * u[1]);
  result[1] = u[1] * (2 * u[0] - 1);

  return result;
}

void Problem_P21(FP_Type h, std::ostream &output)
{
  std_tWrapper f(LotkaVolterra, 2);
  FP_Type t0 = 0.0;
  FP_Type t1 = 10.0;

  VectorD2 u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the explicit Euler method
  Euler Method1(f, t0, u0);
  Method1.iterate_up_to(t1, h);
  Method1.print(output);
}

void Problem_P22(FP_Type h, std::ostream &output)
{
  std_tWrapper f(LotkaVolterra, 2);
  FP_Type t0 = 0.0;
  FP_Type t1 = 10.0;

  VectorD2 u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the classic Runge-Kutta method
  ERK<ERK_04> Method2(f, t0, u0);
  //ERK_Test_O4 Method2(rhs, t0, u0);
  Method2.iterate_up_to(t1, h);
  Method2.print(output);
}

void Test_Sheet2()
{
  std::cout << "Plots for Lotka-Volterra equation (P21), h = 0.1" << std::endl;

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
