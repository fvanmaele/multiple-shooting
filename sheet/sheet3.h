#ifndef SHEET3_H
#define SHEET3_H

#include "../ivp/runge_kutta.h"

void Problem_P32(std::ostream &output1, std::ostream &output2,
                 std::ostream &output3)
{
  RHS_P21 rhs(5, 2, 2, 1);
  FP_Type t0 = 0.0;
  FP_Type t1 = 10.0;

  dealii::Vector<FP_Type> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the Dormand-Prince 45 method
  DOPRI Tableau;
  dealii::LAPACKFullMatrix<FP_Type> A = MFA<49>(7, 7, Tableau.A);

  dealii::Vector<FP_Type> b1(Tableau.b1.begin(), Tableau.b1.end());
  dealii::Vector<FP_Type> b2(Tableau.b2.begin(), Tableau.b2.end());
  dealii::Vector<FP_Type> c(Tableau.c.begin(), Tableau.c.end());

  // Equidistant method
  ERK Equidistant(rhs, t0, u0, A, b1, c);
  Equidistant.iterate_up_to(t1, 1e-3);
  Equidistant.print(output1);
  std::cout << "Amount of steps: " << Equidistant.n_steps() << std::endl;

  // The Dormand-Prince method has order of convergence 5.
  FP_Type EOC1 = eoc_1step(Equidistant, t1, 1e-3);
  std::cout << "EOC: (Dormand-Prince, h = " << 1e-3 << ") "
            << EOC1 << std::endl;

  // Adaptive method of order 5(4)
  ERK Adaptive(rhs, t0, u0, A, b1, b2, c);
  Adaptive.iterate_with_ssc(t1, 1e-1, 1e-4, 5);
  Adaptive.print(output2);

  size_t steps = Adaptive.n_steps();
  size_t misfires = Adaptive.n_misfires();
  std::cout << "Amount of steps: "
            << steps << std::endl
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
