#include <iostream>
#include <string>
#include <getopt.h>
#include <thread>

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
    y[1] = -1./50;

    return y;
  }
};

template <typename Vector>
Vector RHS_ThomasFermi(typename Vector::value_type t, const Vector &u)
{
  Vector y(2);
  y[0] = t * u[1];
  y[1] = 4 * std::pow(u[0], 1.5);

  return y;
}

void Solve_ThomasFermi(int n, bool adaptive)
{
  std_tWrapper f(RHS_ThomasFermi<VectorD2>, 2);
  FP_Type a = 0;
  FP_Type b = 5;

  // Linear separated BVP
  std::vector<FP_Type> A = {1, 0, 0, 0};
  std::vector<FP_Type> B = {0, 0, 1, 0};
  std::vector<FP_Type> c = {1, 0};

  // Define approximate solution for BVP
  CurveTF* eta = new CurveTF;
  assert((*eta)(0)[0] == c[0]);
  assert((*eta)(5)[0] == c[1]);

  // Compute starting trajectory
  std::vector<FP_Type> x;
  if (adaptive)
    // The TOL chosen here should match the initial step width
    // chosen for adaptive methods (1e-3)
    x = trajectory(a, b, f, eta, 1e-3, 2, false);
  else
    x = linspace(a, b, n+1);

  std::cout << "Amount of intervals: "
            << x.size()-1 << std::endl;

  // Starting values s(0)_1 ... s(0)_n
  std::vector<VectorD2> s0;
  for (auto &c : x)
    s0.push_back((*eta)(c));

  // Plot trajectory
  std::ofstream output_file;
  GnuPlot P("TF_trajectory.dat", output_file);

  for (size_t i = 0; i < x.size(); i++)
    output_file << x.at(i) << "\t" << s0.at(i);

  P.plot_with_lines(2, "linespoints");

  // Compute the n IVP's y(t; t_k, s_k)
  SF_External S(f);
  std::vector<VectorD2> y;
  std::vector<VectorD2> s = s0;

  for (size_t i = 0; i < x.size()-1; i++)
    {
      std::cout << x.at(i) << "\t" << x.at(i+1) << "\t" << s.at(i);
      VectorD2 y_stage = S.solve_y(x.at(i), x.at(i+1), s.at(i));
    }

  // ...
}

void usage()
{
  std::cerr << "usage: multiple-shooting [--run-tests] [--intervals <n>] [--adaptive]" << std::endl;
  std::exit(1);
}

int main(int argc, char* argv[])
{
  bool waiting_for_int = false;

  for (int i = 1; i < argc; i++)
    {
      std::string option = argv[i];

      if (waiting_for_int)
        {
          n_intervals = std::stoi(option);
          waiting_for_int = false;
        }
      else if (option == "--run-tests")
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
      else if (option == "--adaptive")
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

  // Linear separated BVP (::LinearBVP)
  Solve_ThomasFermi(n_intervals, adaptive_intervals);

  return 0;
}
