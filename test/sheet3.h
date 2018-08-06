#ifndef SHEET3_H
#define SHEET3_H

#include "../test/sheet2.h"
#include "../ivp/runge_kutta.h"

void Problem_P32(std::ostream &output1, std::ostream &output2,
                 std::ostream &output3)
{
  std_tWrapper f(LotkaVolterra, 2);
  FP_Type t0 = 0.0;
  FP_Type t1 = 30.0;

  VectorD2 u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the Dormand-Prince 87 method
  ERK<DOPRI87> Equidistant(f, t0, u0);
  Equidistant.iterate_up_to(t1, 1e-3);
  Equidistant.print(output1);
  std::cout << "Amount of steps: (DOPRI, Equidistant) "
            << Equidistant.n_steps() << std::endl;

  // Adaptive method of order 8(7)
  ERK<DOPRI87> Adaptive(f, t0, u0);
  FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());
  Adaptive.iterate_with_ssc(t1, 1e-1, TOL, false);
  Adaptive.print(output2);

  size_t steps = Adaptive.n_steps();
  size_t misfires = Adaptive.n_misfires();
  std::cout << "Amount of steps: (DOPRI, Adaptive) "
            << Adaptive.n_steps() << std::endl
            << "Amount of misfires: "
            << misfires << std::endl
            << "Amount of accepted steps: "
            << steps - misfires << std::endl;

  // Plot the adaptive step size on the interval I
  Adaptive.print_step_size(output3);
}

void Test_Sheet3()
{
  std::ofstream output_file1, output_file2, output_file3;
  output_file1.open("volterra_dopri_eq.dat");
  assert(output_file1.is_open());
  output_file2.open("volterra_dopri_ad.dat");
  assert(output_file2.is_open());
  output_file3.open("volterra_dopri_steps.dat");
  assert(output_file3.is_open());

  Problem_P32(output_file1, output_file2, output_file3);
  std::system("gnuplot -p -e \"plot 'volterra_dopri_eq.dat' using 1:2 with lines, "
              "'volterra_dopri_eq.dat' using 1:3 with lines\"");
  std::system("gnuplot -p -e \"plot 'volterra_dopri_ad.dat' using 1:2 with lines, "
              "'volterra_dopri_ad.dat' using 1:3 with lines\"");
  std::system("gnuplot -p -e \"plot 'volterra_dopri_steps.dat' using 1:2 with lines, "
              "'volterra_dopri_steps.dat' using 1:3 with lines\"");
}

#endif // SHEET3_H
