#ifndef SHEET2_H
#define SHEET2_H

void Problem_P21(NumberType h, std::ostream &output)
{
  RHS_P21 rhs(5, 2, 2, 1);
  NumberType t0 = 0.0;
  NumberType t1 = 10.0;

  dealii::Vector<NumberType> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the explicit Euler method
  Explicit_Euler Method1(rhs, t0, u0, h);
  Method1.iterate_up_to(t1);
  Method1.print(output);
}

void Problem_P22(NumberType h, std::ostream &output)
{
  RHS_P21 rhs(5, 2, 2, 1);
  NumberType t0 = 0.0;
  NumberType t1 = 10.0;

  dealii::Vector<NumberType> u0(2);
  u0[0] = 1.0;
  u0[1] = 1.0;

  // Solve the equations with the classic Runge-Kutta method
  ERK_04 Tableau;
  dealii::FullMatrix<NumberType> A(4, 4, Tableau.entries);

  // Use Vector(const InputIterator...) ?
  dealii::Vector<NumberType> b = std_to_Vector(Tableau.b);
  dealii::Vector<NumberType> c = std_to_Vector(Tableau.c);

  ERK Method2(rhs, t0, u0, A, b, c, h);
  Method2.iterate_up_to(t1);
  Method2.print(output);
}

void run_sheet2()
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
