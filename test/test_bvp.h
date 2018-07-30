#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>

#include "../algo/newton.h"
#include "../base/types.h"
#include "../bvp/shooting.h"

// RHS for the initial value problem:
//    w'' = 1.5 * w^2
// resp. the system
//    (u1, u2)' = (u2, 1.5 * u1^2)
class RHS_BVP1 : public RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == 2);
    dealii::Vector<FP_Type> result(2);

    result[0] = u[1];
    result[1] = 1.5 * std::pow(u[0], 2);

    return result;
  }
};

dealii::Vector<FP_Type>
perform_step(DivFunctor &F, dealii::Vector<FP_Type> s)
{
  Newton N(F, s.size());
  size_t steps = 0;

  for (size_t k = 0; k < 50; k++)
    {
      // 2. Compute approximate Jacobian of F in s_k
      // 3. Perform step of (quasi-)Newton method
      s = N.step(F.jacobian(s), s, true);
      steps++;

      if (N.stopping_criterion(1e-8))
        {
          std::cout << "Solution of F: (" << steps << " steps) " << s;
          return s;
        }
    }
}

// Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
void Test_ED()
{
  RHS_BVP1 f;
  FP_Type a = 0.0;
  FP_Type b = 1.0;

  // linear separated boundary value problem
  std::vector<FP_Type> _A = { 1, 0,
                              0, 0 };
  std::vector<FP_Type> _B = { 0, 0,
                              1, 0 };
  dealii::FullMatrix A(2, 2, _A.data());
  dealii::FullMatrix B(2, 2, _B.data());

  dealii::Vector<FP_Type> c(2);
  c[0] = 4.;
  c[1] = 1.;

  // F used in Newton iteration (root indicates a solution
  // to the boundary value problem)
  ShootingFunction_ED F(f, a, b, A, B, c);

  // starting value
  dealii::Vector<FP_Type> s(2);
  s[0] = 4;
  s[1] = -100;

  // graph of F_1(s_1)
  std::ofstream output_file;
  output_file.open("bvp_sval.dat");
  assert(output_file.is_open());

  // s = [4, -100] .. [4, 0], step width 1
  for (size_t i = 0; i < 100; i++)
    {
      s[1] += 1;
      dealii::Vector<FP_Type> diff = F.value(s);

      output_file << s[1] << "\t" << diff[1] << std::endl;
    }
  std::system("gnuplot -p -e \"plot 'bvp_sval.dat' using 1:2 with lines, "
              "'bvp_sval.dat' using 1:3 with lines\"");

  // Newton method
  s[0] = 4;
  s[1] = -30;
  dealii::Vector<FP_Type> u = perform_step(F, s);

  s[1] = -1;
  dealii::Vector<FP_Type> v = perform_step(F, s);
}

#endif // BVP_EXTERNAL_H
