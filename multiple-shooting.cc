#include <iostream>
#include <string>
#include <getopt.h>
#include <thread>

#include <deal.II/lac/block_vector.h>

#include "base/types.h"
#include "bvp/boundary.h"
#include "bvp/methods.h"
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

void Solve_ThomasFermi(int n_int, bool adaptive)
{
  std_tWrapper f(RHS_ThomasFermi<VectorD2>, 2);
  const FP_Type a = 0;
  const FP_Type b = 5;

  // Linear separated BVP
  const MatrixD2 A = init_matrix(2, 2, {1, 0, 0, 0});
  const MatrixD2 B = init_matrix(2, 2, {0, 0, 1, 0});
  const VectorD2 c = init_vector(2, {1, 0});
  BC_Linear r(A, B, c);

  // Define approximate solution for BVP
  CurveTF* eta = new CurveTF;
  assert((*eta)(0)[0] == c[0]);
  assert((*eta)(5)[0] == c[1]);

  // Subdivide interval a = t0 < ... < t1 = b
  std::vector<FP_Type> t;

  if (adaptive)
    // The TOL chosen here should match the initial step width
    // chosen for adaptive methods (1e-3)
    t = trajectory(a, b, f, eta, 1e-3, 2, false);
  else
    t = linspace(a, b, n_int+1);

  assert(t.front() == a);
  assert(t.back()  == b);

  const size_t n = 2;
  const size_t m = t.size();
  std::cout << "Amount of intervals: " << m-1 << std::endl;

  // Starting values s(0)_1 ... s(0)_n
  std::vector<VectorD2> s0;

  for (auto &c : t)
    s0.push_back((*eta)(c));

  // Plot trajectory
  std::ofstream output_file;
  GnuPlot Dat1("TF_trajectory.dat", output_file);

  for (size_t i = 0; i < t.size(); i++)
    output_file << t.at(i) << "\t" << s0.at(i);
  Dat1.plot_with_lines(2, "linespoints");

  // Begin multiple shooting method
  // ...
}

void usage()
{
  std::cerr << "usage: multiple-shooting [--run-tests] [--intervals <n>] [--adaptive]"
            << std::endl;
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
