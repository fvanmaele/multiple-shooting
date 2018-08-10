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
  // Time interval
  const FP_Type a = 0;
  const FP_Type b = 5;

  // RHS of ODE
  const size_t dim = 2;
  FAD_tWrapper f(RHS_ThomasFermi<VectorAD>, dim);

  // Boundary conditions (linear separated BVP)
  const MatrixD2 A = init_matrix(2, 2, {1, 0, 0, 0});
  const MatrixD2 B = init_matrix(2, 2, {0, 0, 1, 0});
  const VectorD2 c = init_vector(2, {1, 0});
  BC_Linear r(A, B, c);

  // Approximate solution for BVP
  CurveTF* eta = new CurveTF;
  assert((*eta)(0)[0] == c[0]);
  assert((*eta)(5)[0] == c[1]);

  // 1) Subdivide interval a = t0 < ... < t1 = b
  std::vector<FP_Type> t;

  if (adaptive)
    // The TOL chosen here should match the initial step width
    // chosen for adaptive methods (1e-3)
    t = trajectory(a, b, f, eta, 1e-3, 2, false);
  else
    t = linspace(a, b, n_int+1);

  if (!(t.front() == a && t.back() == b))
    throw std::domain_error("subdivision does not match interval boundaries");

  const size_t m = t.size();
  std::cout << "Amount of intervals: " << m-1 << std::endl;

  // 2) Starting values s(0)_1 ... s(0)_n
  VectorD2 s0(m*dim);
  for (size_t i = 0; i < m; i++)
    { // Extract i-th block of s
      VectorD2 s_i = (*eta)(t.at(i));

      for (size_t k = 0; k < dim; k++)
        {
          s0[k + i*dim] = s_i[k];
        }
    }

  // a) Plot trajectory
  std::ofstream output_file;
  GnuPlot Dat1("TF_trajectory.dat", output_file);

  for (size_t i = 0; i < t.size(); i++)
    {
      VectorD2 s_i = (*eta)(t.at(i));
      output_file << t.at(i) << "\t" << s_i;
    }
  Dat1.plot_with_lines(2, "linespoints");

  // 3) Begin multiple shooting method
  MultipleShooting<SF_External> F(f, dim, t, r);
  Newton N(F, m * dim);
  VectorD2 sol = N.iterate(s0);

  // 4) Plot graph of solution

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
