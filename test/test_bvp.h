#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

#include "../base/types.h"
#include "../base/gnuplot.h"
#include "../lac/lac_types.h"
#include "../bvp/linear.h"
#include "../bvp/trajectory.h"

namespace Test
{
  // RHS for the initial value problem:
  //    w'' = 1.5 * w^2   resp.   (u1, u2)' = (u2, 1.5 * u1^2)
  template <typename Vector>
  Vector RHS_Stoer(typename Vector::value_type, const Vector &u)
  {
    Vector result(2);
    result[0] = u[1];
    result[1] = 1.5 * std::pow(u[0], 2);

    return result;
  }

  template <typename Vector>
  Vector RHS_Troesch(typename Vector::value_type, const Vector &u)
  {
    Vector result(2);
    result[0] = u[1];
    result[1] = 5 * std::sinh(5 * u[0]);

    return result;
  }

  class CurveStoer : public Curve
  {
  public:
    VectorD2 operator()(FP_Type t)
    {
      VectorD2 y(2);
      y[0] = -3*t + 4;
      y[1] = -3;

      return y;
    }
  };

  template <typename DiffMethod>
  void Stoer_Graph(LinearBVP<DiffMethod> &F, std::string filename)
  {
    VectorD2 s(2);
    s[0] = 4;
    s[1] = -100;

    std::vector<VectorD2> range;
    for (size_t i = 0; i < 100; i++)
      {
        s[1] += 1;
        range.emplace_back(s);
      }

    std::ofstream output_file;
    GnuPlot Dat1(filename, output_file);

    for (auto &s : range)
      output_file << s[1] << "\t" << F(s)[1] << std::endl;
    Dat1.plot_with_lines(1);
  }

  void Stoer_Mult(TimeFunctor &rhs, FP_Type a, FP_Type b)
  {
    CurveStoer* eta = new CurveStoer;
    assert((*eta)(0)[0] == 4);
    assert((*eta)(1)[0] == 1);

    std::vector<FP_Type> subint = trajectory(a, b, rhs, eta, 1e-3, 2, false);
    std::cout << "Amount of intervals: " << subint.size()-1 << std::endl;
    std::cout << subint;

    // Plot trajectory
    std::ofstream output_file;
    GnuPlot Dat1("Stoer_trajectory.dat", output_file);

    for (auto &c : subint)
      output_file << c << "\t" << (*eta)(c);
    Dat1.plot_with_lines(2, "linespoints");
  }

  void Stoer()
  {
    // Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
    std_tWrapper rhs(RHS_Stoer<VectorD2>, 2);
    FAD_tWrapper rhs_ad(RHS_Stoer<VectorAD>, 2);

    FP_Type a = 0.0;
    FP_Type b = 1.0;

    // Linear separated BVP
    std::vector<FP_Type> A = {1, 0, 0, 0};
    std::vector<FP_Type> B = {0, 0, 1, 0};
    std::vector<FP_Type> c = {4, 1};

    // Initial values s0_1, s0_2 for solutions s1, s2
    std::vector<VectorD2> start;
    VectorD2 s(2);

    s[0] = 4;
    s[1] = -99;
    start.emplace_back(s);

    s[0] = 4;
    s[1] = -1;
    start.emplace_back(s);

    // Solve BVP (external differentation)
    std::cout << "Single shooting (Stoer, ext. diff.)" << std::endl;
    LinearBVP<SF_External> F(rhs, a, b, A, B, c);
    Newton N(F, 2);

    for (auto &s : start)
      N.iterate(s);

    // Solve BVP (automatic differentation)
    std::cout << "Single shooting (Stoer, aut. diff.)" << std::endl;
    LinearBVP<SF_Automatic> F_ad(rhs_ad, a, b, A, B, c);
    Newton N_ad(F_ad, 2);

    for (auto &s : start)
      N_ad.iterate(s);

    // Create plot for F(s)
    Stoer_Graph<SF_External>(F, "bvp_stoer_ed.dat");
    Stoer_Graph<SF_Automatic>(F_ad, "bvp_stoer_ad.dat");

    // Find subintervals for multiple shooting method
    Stoer_Mult(rhs, a, b);
  }

  class CurveTroesch : public Curve
  {
  public:
    VectorD2 operator()(FP_Type t)
    {
      VectorD2 y(2);
      y[0] = t;
      y[1] = 1;

      return y;
    }
  };

  void Troesch_Mult(TimeFunctor &rhs, FP_Type a, FP_Type b)
  {
    // Interval subdivision for multiple shooting method
    CurveTroesch* eta = new CurveTroesch;
    assert((*eta)(0)[0] == 0);
    assert((*eta)(1)[0] == 1);

    std::vector<FP_Type> subint = trajectory(a, b, rhs, eta, 1e-3, 2, false);
    std::cout << "Amount of intervals: " << subint.size()-1 << std::endl;
    std::cout << subint;

    // Plot trajectory
    std::ofstream output_file;
    GnuPlot Dat1("Troesch_trajectory.dat", output_file);

    for (auto &c : subint)
      output_file << c << "\t" << (*eta)(c);
    Dat1.plot_with_lines(2, "linespoints");
  }

  void Troesch()
  {
    // Troesch BVP, see Num. Math. 2, 7.4.3.8
    std_tWrapper rhs(RHS_Troesch<VectorD2>, 2);
    FAD_tWrapper rhs_ad(RHS_Troesch<VectorAD>, 2);

    FP_Type a = 0.0;
    FP_Type b = 1.0;

    // Linear separated BVP
    std::vector<FP_Type> A = {1, 0, 0, 0};
    std::vector<FP_Type> B = {0, 0, 1, 0};
    std::vector<FP_Type> c = {0, 1};

    // Solve BVP (external differentation)
    std::cout << "Single shooting (Troesch, ext. diff.)" << std::endl;
    LinearBVP<SF_External> F(rhs, a, b, A, B, c);
    Newton N(F, 2);

    // Solve BVP (automatic differentation)
    std::cout << "Single shooting (Troesch, aut. diff.)" << std::endl;
    LinearBVP<SF_Automatic> F_ad(rhs_ad, a, b, A, B, c);
    Newton N_ad(F_ad, 2);

    // Left solution
    VectorD2 s(2);
    s[0] = 0;
    s[1] = 0.05;

    N.iterate(s);
    N_ad.iterate(s);

    // Find subintervals for multiple shooting method
    Troesch_Mult(rhs, a, b);
  }
}

#endif // BVP_EXTERNAL_H
