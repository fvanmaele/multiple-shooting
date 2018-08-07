#ifndef SHEET3_H
#define SHEET3_H

#include "../test/sheet2.h"
#include "../ivp/runge_kutta.h"

namespace Test
{
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
    ERK<DOPRI> Equidistant(f, t0, u0);
    Equidistant.iterate_up_to(t1, 1e-3);
    Equidistant.print(output1);
    std::cout << "Amount of steps: (DOPRI54, Equidistant) "
              << Equidistant.n_steps() << std::endl;

    // Adaptive method of order 5(4)
    ERK<DOPRI> Adaptive54(f, t0, u0);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());
    Adaptive54.iterate_with_ssc(t1, 1e-1, TOL, false);

    // Adaptive method of order 8(7)
    ERK<DOPRI87> Adaptive87(f, t0, u0);
    Adaptive87.iterate_with_ssc(t1, 1e-1, TOL, false);
    Adaptive87.print(output2);

    std::cout << "Amount of steps: (DOPRI54, Adaptive) "
              << Adaptive54.n_steps() << std::endl
              << "Amount of misfires: "
              << Adaptive54.n_misfires() << std::endl
              << "Amount of accepted steps: "
              << Adaptive54.n_steps() - Adaptive54.n_misfires() << std::endl
              << "Amount of steps: (DOPRI87, Adaptive) "
              << Adaptive87.n_steps() << std::endl
              << "Amount of misfires: "
              << Adaptive87.n_misfires() << std::endl
              << "Amount of accepted steps: "
              << Adaptive87.n_steps() - Adaptive87.n_misfires() << std::endl;

    // Plot the adaptive step size on the interval I
    Adaptive87.print_step_size(output3);
  }

  void Sheet3()
  {
    std::ofstream output_file1, output_file2, output_file3;
    output_file1.open("volterra_dopri54_eq.dat");
    assert(output_file1.is_open());
    output_file2.open("volterra_dopri87_ad.dat");
    assert(output_file2.is_open());
    output_file3.open("volterra_dopri87_steps.dat");
    assert(output_file3.is_open());

    Problem_P32(output_file1, output_file2, output_file3);
    std::system("gnuplot -p -e \"plot 'volterra_dopri54_eq.dat' using 1:2 with lines, "
                "'volterra_dopri54_eq.dat' using 1:3 with lines\"");
    std::system("gnuplot -p -e \"plot 'volterra_dopri87_ad.dat' using 1:2 with lines, "
                "'volterra_dopri87_ad.dat' using 1:3 with lines\"");
    std::system("gnuplot -p -e \"plot 'volterra_dopri87_steps.dat' using 1:2 with lines, "
                "'volterra_dopri87_steps.dat' using 1:3 with lines\"");
  }
}

#endif // SHEET3_H
