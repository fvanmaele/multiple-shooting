#include <chrono>
#include <iostream>
#include <string>

#include "base/types.h"
#include "bvp/boundary.h"
#include "bvp/methods.h"
#include "bvp/shooting.h"
#include "bvp/trajectory.h"
#include "ivp/runge_kutta.h"
#include "test/run_tests.h"

static int n_intervals = 20;
static bool run_tests = false;
static bool use_broyden_method = false;

class CurveTF : public Curve
{
public:
  using Curve::Curve;

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

void Solve_ThomasFermi(int n_int)
{
  auto start0 = std::chrono::system_clock::now();
  // -------------------------------------------
  // 0) Linear separated BVP
  const FP_Type a = 0;
  const FP_Type b = 5;

  // Boundary conditions
  const MatrixD2 A = init_matrix(2, 2, {1, 0, 0, 0});
  const MatrixD2 B = init_matrix(2, 2, {0, 0, 1, 0});
  const VectorD2 c = init_vector(2, {1, 0});
  BC_Linear r(A, B, c);

  // Approximate solution (starting trajectory)
  CurveTF eta(2);
  assert(eta(0)[0] == c[0]);
  assert(eta(5)[0] == c[1]);

  // RHS of ODE
  const size_t dim = eta.n_dim();
  FAD_tWrapper f(RHS_ThomasFermi<VectorAD>, dim);

  
  auto start1 = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds01 = start1 - start0;
  std::cout << "STEP 0: elapsed time: " << elapsed_seconds01.count() << "s" << std::endl;

  // -------------------------------------------
  // 1) Subdivide interval a = t0 < ... < t1 = b
  std::vector<FP_Type> t;

  // a) Find initial subdivision based on starting trajectory
  t = trajectory(a, b, f, &eta, 1.1, false, 1e-2);

  if (t.size() < (size_t)(n_int+1))
    throw std::invalid_argument("insufficient time points available");

  // b) Interpolate to given interval amount
  t = interpolate_points(t, n_int+1);
  std::cout << "Amount of intervals: " << t.size()-1 << std::endl;
  std::cout << t;

  // c) Check boundaries
  if (!(t.front() == a && t.back() == b))
    throw std::domain_error("subdivision does not match interval boundaries");

  // d) Plot trajectory
  std::ofstream output_file1;
  GnuPlot Dat1("TF_trajectory.dat", output_file1);

  for (auto &c : t)
    output_file1 << c << "\t" << eta(c);


  auto start2 = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds12 = start2 - start1;
  std::cout << "STEP 1: elapsed time: " << elapsed_seconds12.count() << "s" << std::endl;
  
  // -------------------------------------------
  // 2) Starting values s(0)_1 ... s(0)_n
  const size_t m = t.size();
  VectorD2 s0(m * dim);

  for (size_t i = 0; i < m; i++)
    { // Extract i-th block of s
      VectorD2 s_i = eta(t.at(i));

      for (size_t k = 0; k < dim; k++)
        s0[k + i*dim] = s_i[k];
    }

  
  auto start3 = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds23 = start3 - start2;
  std::cout << "STEP 2: elapsed time: " << elapsed_seconds23.count() << "s" << std::endl;
  
  // -------------------------------------------
  // 3) Begin multiple shooting method
  SF_Automatic<KARP> M(f, true, 1e-2, 1e-4);
  MultipleShooting F(M, t, r);
  Newton N(F, m * dim);

  VectorD2 sol;
  if (use_broyden_method)
    sol = N.iterate_broyden(s0);
  else
    sol = N.iterate(s0);


  auto start4 = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds34 = start4 - start3;
  std::cout << "STEP 3: elapsed time: " << elapsed_seconds34.count() << "s" << std::endl;
  
  // -------------------------------------------
  // 4) Plot graph of solution
  std::ofstream output_file2;
  GnuPlot Dat2("TF_solution.dat", output_file2);

  for (size_t i = 0; i < t.size()-1; i++)
    {
      VectorD2 s_i(dim);

      for (size_t k = 0; k < dim; k++)
        s_i[k] = sol[k + i*dim];

      // Solve u' = f(t_{i+1}; t_i, s_i)
      ERK<KARP> M_i(f, t[i], s_i);
      M_i.iterate_with_ssc(t[i+1], 1e-2, 1e-7);
      M_i.print(output_file2);
    }


  auto start5 = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds45 = start5 - start4;
  std::cout << "STEP 4: elapsed time: " << elapsed_seconds45.count() << "s" << std::endl;

  // -------------------------------------------
  // 5) Plot graphs
  Dat1.plot_with_lines(2, "linespoints");
  Dat2.plot_with_lines(2, "linespoints");
  Dat2.plot_with_lines(2, "lines", true);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<FP_Type> elapsed_seconds_end = end - start5;
  std::cout << "GNUPLOT: elapsed time: " << elapsed_seconds_end.count() << "s" << std::endl;

  std::chrono::duration<FP_Type> elapsed_seconds_total = end - start0;
  std::cout << "Total elapsed time: " << elapsed_seconds_total.count() << "s" << std::endl;
}

void usage()
{
  std::cerr << "usage: multiple-shooting [--run-tests] [--intervals <n>] [--broyden]"
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
      else if (option == "--broyden")
        {
          use_broyden_method = true;
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

  if (n_intervals < 0)
    {
      std::cerr << "interval number must be positive" << std::endl;
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
  Solve_ThomasFermi(n_intervals);

  return 0;
}
