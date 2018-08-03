#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../bvp/linear.h"

// RHS for the initial value problem:
//    w'' = 1.5 * w^2
//
// resp. the system
//    (u1, u2)' = (u2, 1.5 * u1^2)
template <typename Vector>
Vector Stoer(typename Vector::value_type t, const Vector &u)
{
  Vector result(2);
  result[0] = u[1];
  result[1] = 1.5 * std::pow(u[0], 2);

  return result;
}

template <typename Vector>
Vector Troesch(typename Vector::value_type t, const Vector &u)
{
  Vector result(2);
  result[0] = u[1];
  result[1] = 5 * std::sinh(5 * u[0]);

  return result;
}

// Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
void Test_Stoer()
{
  std_tWrapper f(Stoer<VectorD2>, 2);
  FAD_tWrapper f_auto(Stoer<VectorAD>, 2);
  FP_Type a = 0.0;
  FP_Type b = 1.0;

  dealii::Vector<FP_Type> c(2);
  c[0] = 4.;
  c[1] = 1.;

  std::vector<dealii::Vector<FP_Type> > start;
  dealii::Vector<FP_Type> s(2);

  // Left solution
  s[0] = 4;
  s[1] = -99;
  start.emplace_back(s);

  // Right solution
  s[0] = 4;
  s[1] = -1;
  start.emplace_back(s);

  std::cout << "Single shooting (Stoer) - External differentation" << std::endl;
  SimpleBVP<SF_External> BVP(f, a, b, c);
  BVP.single_shooting(start);

  std::cout << "Single shooting (Stoer) - Automatic differentation" << std::endl;
  SimpleBVP<SF_Automatic> BVP_auto(f_auto, a, b, c);
  BVP_auto.single_shooting(start);

  // Create plot
  std::ofstream output_file1, output_file2;
  output_file1.open("bvp_sval_ed.dat");
  assert(output_file1.is_open());
  output_file2.open("bvp_sval_ad.dat");
  assert(output_file2.is_open());

  std::vector<dealii::Vector<FP_Type> > range;
  s[0] = 4;
  s[1] = -100;

  for (size_t i = 0; i < 100; i++)
    {
      s[1] += 1;
      range.emplace_back(s);
    }
  BVP.shooting_graph(2, range, output_file1);
  BVP_auto.shooting_graph(2, range, output_file2);

  std::system("gnuplot -p -e \"plot 'bvp_sval_ed.dat' using 2:4 with lines\"");
  std::system("gnuplot -p -e \"plot 'bvp_sval_ad.dat' using 2:4 with lines\"");
}

// Troesch BVP, see Num. Math. 2, 7.4.3.8
void Test_Troesch()
{
  std_tWrapper f(Troesch<VectorD2>, 2);
  FAD_tWrapper f_auto(Troesch<VectorAD>, 2);
  FP_Type a = 0.0;
  FP_Type b = 1.0;

  dealii::Vector<FP_Type> c(2);
  c[0] = 0.;
  c[1] = 1.;

  std::vector<dealii::Vector<FP_Type> > start;
  dealii::Vector<FP_Type> s(2);

  // Left solution
  s[0] = 0;
  s[1] = 0.05;
  start.emplace_back(s);

  std::cout << "Single shooting (Troesch) - External differentation" << std::endl;
  SimpleBVP<SF_External> BVP(f, a, b, c);
  BVP.single_shooting(start);

  std::cout << "Single shooting (Troesch) - Automatic differentation" << std::endl;
  SimpleBVP<SF_Automatic> BVP_auto(f_auto, a, b, c);
  BVP_auto.single_shooting(start);
}

#endif // BVP_EXTERNAL_H
