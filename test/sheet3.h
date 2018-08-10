#ifndef SHEET3_H
#define SHEET3_H

#include "../base/gnuplot.h"
#include "../test/sheet2.h"
#include "../ivp/runge_kutta.h"

namespace Test
{
  void Problem_P32(std::ostream &output1, std::ostream &output2,
                   std::ostream &output3)
  {
    std_tWrapper f(LotkaVolterra, 2);
    FP_Type t0 = 0.0;
    FP_Type t1 = 20.0;

    VectorD2 u0(2);
    u0[0] = 1.0;
    u0[1] = 1.0;

    // Solve the equations with the Dormand-Prince 87 method
    ERK<DOPRI54> Equidistant(f, t0, u0);
    Equidistant.iterate_up_to(t1, 1e-3);
    Equidistant.print(output1);
    std::cout << "Amount of steps: (DOPRI54, Equidistant) "
              << Equidistant.n_steps() << std::endl;

    // Adaptive method of order 5(4)
    ERK<DOPRI54> Adaptive54(f, t0, u0);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());
    Adaptive54.iterate_with_ssc(t1, 1e-1, TOL, false);
    Adaptive54.print(output2);

    std::cout << "Amount of steps: (DOPRI54, Adaptive) "
              << Adaptive54.n_steps() << std::endl
              << "Amount of misfires: "
              << Adaptive54.n_misfires() << std::endl
              << "Amount of accepted steps: "
              << Adaptive54.n_steps() - Adaptive54.n_misfires() << std::endl;

    // Plot the adaptive step size on the interval I
    Adaptive54.print_step_size(output3);
  }

  void Sheet3()
  {
    std::ofstream output_file1, output_file2, output_file3;
    GnuPlot Dat1("volterra_dopri54_eq.dat", output_file1);
    GnuPlot Dat2("volterra_dopri54_ad.dat", output_file2);
    GnuPlot Dat3("volterra_dopri54_steps.dat", output_file3);

    Problem_P32(output_file1, output_file2, output_file3);
    Dat1.plot_with_lines(2);
    Dat2.plot_with_lines(2);
    Dat3.plot_with_lines(2);
  }
}

#endif // SHEET3_H
