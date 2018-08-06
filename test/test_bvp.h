#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../bvp/linear.h"

// RHS for the initial value problem:
//    w'' = 1.5 * w^2
//
// resp. the system
//    (u1, u2)' = (u2, 1.5 * u1^2)
template <typename Vector>
Vector Stoer(typename Vector::value_type, const Vector &u)
{
  Vector result(2);
  result[0] = u[1];
  result[1] = 1.5 * std::pow(u[0], 2);

  return result;
}

template <typename Vector>
Vector Troesch(typename Vector::value_type, const Vector &u)
{
  Vector result(2);
  result[0] = u[1];
  result[1] = 5 * std::sinh(5 * u[0]);

  return result;
}

template <typename DiffMethod>
void Stoer_Graph(SimpleBVP<DiffMethod> &BVP, std::string filename)
{
  // Create plot
  std::ofstream output_file1;
  output_file1.open(filename.c_str());
  assert(output_file1.is_open());

  std::vector<dealii::Vector<FP_Type> > range;
  dealii::Vector<FP_Type> s(2);
  s[0] = 4;
  s[1] = -100;

  for (size_t i = 0; i < 100; i++)
    {
      s[1] += 1;
      range.emplace_back(s);
    }

  BVP.shooting_graph(2, range, output_file1);
  std::string cmd = "gnuplot -p -e \"plot '" + filename + "' using 2:4 with lines\"";
  std::system(cmd.c_str());
}

// Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
void Test_Stoer()
{
  std_tWrapper f(Stoer<VectorD2>, 2);
  FAD_tWrapper f_ad(Stoer<VectorAD>, 2);
  FP_Type a = 0.0;
  FP_Type b = 1.0;

  VectorD2 c(2);
  c[0] = 4.;
  c[1] = 1.;

  std::vector<VectorD2> start;
  VectorD2 s(2);

  // Left solution
  s[0] = 4;
  s[1] = -99;
  start.emplace_back(s);

  // Right solution
  s[0] = 4;
  s[1] = -1;
  start.emplace_back(s);

  std::cout << "Single shooting (Stoer, ext. diff.)" << std::endl;
  SimpleBVP<SF_External> BVP(f, a, b, c);
  BVP.single_shooting(start);

  std::cout << "Single shooting (Stoer, aut. diff.)" << std::endl;
  SimpleBVP<SF_Automatic> BVP_AD(f_ad, a, b, c);
  BVP_AD.single_shooting(start);

  Stoer_Graph<SF_External>(BVP, "bvp_stoer_ed.dat");
  Stoer_Graph<SF_Automatic>(BVP_AD, "bvp_stoer_ad.dat");
}

// Troesch BVP, see Num. Math. 2, 7.4.3.8
void Test_Troesch()
{
  std_tWrapper f(Troesch<VectorD2>, 2);
  FAD_tWrapper f_ad(Troesch<VectorAD>, 2);
  FP_Type a = 0.0;
  FP_Type b = 1.0;

  VectorD2 c(2);
  c[0] = 0.;
  c[1] = 1.;

  std::vector<VectorD2> start;
  VectorD2 s(2);

  // Left solution
  s[0] = 0;
  s[1] = 0.05;
  start.emplace_back(s);

  std::cout << "Single shooting (Troesch, ext. diff.)" << std::endl;
  SimpleBVP<SF_External> BVP(f, a, b, c);
  BVP.single_shooting(start);

  std::cout << "Single shooting (Troesch, aut. diff.)" << std::endl;
  SimpleBVP<SF_Automatic> BVP_AD(f_ad, a, b, c);
  BVP.single_shooting(start);
}

#endif // BVP_EXTERNAL_H
