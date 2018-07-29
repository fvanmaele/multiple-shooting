#ifndef SHEET3_H
#define SHEET3_H

#include "../ivp/runge_kutta.h"

void Problem_P32(std::ostream &output1, std::ostream &output2,
                 std::ostream &output3)
{
  RHS_P21 rhs(5, 2, 2, 1);
  FP_Type t0 = 0.0;
  FP_Type t1 = 50.0;

  dealii::Vector<FP_Type> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the Dormand-Prince 45 method
  DOPRI T;
  ButcherTableau<7> Tab(T.A, T.b1, T.b2, T.c);

  // Equidistant method
  ERK Equidistant(rhs, t0, u0, Tab.matrix, Tab.weights, Tab.nodes);
  Equidistant.iterate_up_to(t1, 1e-3);
  Equidistant.print(output1);
  std::cout << "Amount of steps: " << Equidistant.n_steps() << std::endl;

  // Adaptive method of order 5(4)
  ERK Adaptive(rhs, t0, u0, Tab.matrix, Tab.weights, Tab.weights_low, Tab.nodes);
  Adaptive.iterate_with_ssc(t1, 1e-1, 1e-8, 4);
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
