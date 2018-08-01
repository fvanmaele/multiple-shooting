#ifndef SHEET2_H
#define SHEET2_H

#include <cassert>

#include "../algo/convergence.h"
#include "../base/types.h"
#include "../ivp/euler.h"
#include "../ivp/runge_kutta.h"

// This functor models the Lotka-Volterra equations.
class RHS_P21 : public TimeFunctor
{
public:
  // Constructor for the 4 parameters a, b, c, d.
  RHS_P21(FP_Type _a, FP_Type _b, FP_Type _c, FP_Type _d) :
    a(_a), b(_b), c(_c), d(_d)
  {}

  // u represents the vector u = (u1(t), u2(t)).
  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    dealii::Vector<FP_Type> result(2);
    result[0] = u[0] * (a - b * u[1]);
    result[1] = u[1] * (c * u[0] - d);

    return result;
  }

private:
  FP_Type a, b, c, d;
};

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
  ERK<ERK_04> Method2(rhs, t0, u0);
  //ERK_Test_O4 Method2(rhs, t0, u0);
  Method2.iterate_up_to(t1, h);
  Method2.print(output);
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
