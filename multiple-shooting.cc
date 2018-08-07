#include <iostream>
#include <string>
#include <getopt.h>

#include "base/types.h"
#include "bvp/shooting.h"
#include "bvp/trajectory.h"
#include "ivp/runge_kutta.h"
#include "test/run_tests.h"

static int n_intervals = 20;
static bool adaptive_intervals = false;
static bool run_tests = false;

class CurveTF : public Curve
{
public:
  VectorD2 operator()(FP_Type t)
  {
    VectorD2 y(2);
    y[0] = -1./5*t + 1;
    y[1] = -1./500;

    return y;
  }
};

VectorD2 ThomasFermi(FP_Type t, const VectorD2 &u)
{
  VectorD2 y(2);
  y[0] = t * u[1];
  y[1] = 4 * std::pow(u[0], 1.5);

  return y;
}

void Solve_ThomasFermi(int n, bool adaptive)
{
  std_tWrapper f(ThomasFermi, 2);
  FP_Type a = 0;
  FP_Type b = 5;

  VectorD2 c(2);
  c[0] = 1;
  c[1] = 0;

  // Define approximate solution for BVP
  CurveTF* eta = new CurveTF;
  assert((*eta)(0)[0] == 1);
  assert((*eta)(5)[0] == 0);

  // Compute starting trajectory
  std::vector<FP_Type> subint;
  if (adaptive)
    subint = trajectory(a, b, f, eta, 1e-2, 1.5, false);
  else
    subint = linspace(a, b, n+1);

  std::cout << "Amount of intervals: " << subint.size()-1 << std::endl;
  std::cout << subint;

  // Plot trajectory
  std::ofstream output_file;
  GnuPlot P("TF_trajectory.dat", output_file);

  for (auto &c : subint)
    output_file << c << "\t" << (*eta)(c);

  P.plot_with_lines(1, "linespoints");
}

void usage()
{
  std::cerr << "usage: multiple-shooting [--run-tests]" << std::endl;
  std::exit(1);
}

int main(int argc, char* argv[])
{
  bool waiting_for_int(false);

  for (int i = 1; i < argc; i++)
    {
      std::string option = argv[i];

      if (waiting_for_int)
        {
          n_intervals = std::stoi(option);
        }
      if (option == "--run-tests")
        {
          std::cerr << "Running self-tests" << std::endl;
          run_tests = true;
        }
      else if (option == "--help")
        {
          usage();
        }
      else if (option == "--intervals")
        {
          waiting_for_int = true;
        }
      else if (option == "--adaptive-intervals")
        {
          adaptive_intervals = true;
        }
      else
        {
          std::cerr << "unknown option" << std::endl;
          usage();
        }
    }

  if (waiting_for_int)
    {
      std::cerr << "--intervals expects an argument" << std::endl;
      usage();
    }

  if (run_tests)
    {
      Test::exercise_sheet();
      Test::aut_diff();
      Test::linear_bvp();
      Test::var_equation();
    }

  // Linear separated BVP (::SimpleBVP)
  Solve_ThomasFermi(n_intervals, adaptive_intervals);

  return 0;
}
