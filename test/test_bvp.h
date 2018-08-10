#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

#include "../base/types.h"
#include "../base/gnuplot.h"
#include "../lac/lac_types.h"
#include "../bvp/boundary.h"
#include "../bvp/methods.h"
#include "../bvp/shooting.h"
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
    using Curve::Curve;

    VectorD2 operator()(FP_Type t)
    {
      VectorD2 y(2);
      y[0] = -3*t + 4;
      y[1] = -3;

      return y;
    }
  };

  void Stoer_Graph(SingleShooting &F, std::string filename)
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
    CurveStoer* eta = new CurveStoer(2);
    assert((*eta)(0)[0] == 4);
    assert((*eta)(1)[0] == 1);

    std::vector<FP_Type> t = trajectory(a, b, rhs, eta, 2, false, 1e-3);
    std::cout << "Amount of intervals: " << t.size()-1 << std::endl;
    std::cout << t;

    // Plot trajectory
    std::ofstream output_file;
    GnuPlot Dat1("Stoer_trajectory.dat", output_file);

    for (auto &c : t)
      output_file << c << "\t" << (*eta)(c);
    Dat1.plot_with_lines(2, "linespoints");
  }

  void Stoer()
  {
    // Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
    FP_Type a = 0.0;
    FP_Type b = 1.0;

    // Linear separated BVP
    MatrixD2 A = init_matrix(2, 2, {1, 0, 0, 0});
    MatrixD2 B = init_matrix(2, 2, {0, 0, 1, 0});
    VectorD2 c = init_vector(2, {4, 1});
    BC_Linear r(A, B, c);

    // Initial values s0_1, s0_2 for solutions s1, s2
    std::vector<VectorD2> start;
    VectorD2 s(2);

    s[0] = 4;
    s[1] = -99;
    start.emplace_back(s);

    s[0] = 4;
    s[1] = -1;
    start.emplace_back(s);

    // 1. single shooting, external differentation
    std::cout << "Single shooting (Stoer, ext. diff.)" << std::endl;
    std_tWrapper rhs(RHS_Stoer<VectorD2>, 2);
    SF_External<DOPRI54> M(rhs, true, 1e-1, 1e-4);
    SingleShooting F(M, a, b, r);
    Newton N(F, 2);

    for (auto &s : start)
      N.iterate(s);

    // 2. single shooting, automatic differentation
    std::cout << "Single shooting (Stoer, aut. diff.)" << std::endl;
    FAD_tWrapper rhs_ad(RHS_Stoer<VectorAD>, 2);
    SF_Automatic<DOPRI54> M_ad(rhs_ad, true, 1e-1, 1e-4);
    SingleShooting F_ad(M_ad, a, b, r);
    Newton N_ad(F_ad, 2);

    for (auto &s : start)
      N_ad.iterate(s);

    // Create plot for F(s)
    Stoer_Graph(F, "bvp_stoer_ed.dat");
    Stoer_Graph(F_ad, "bvp_stoer_ad.dat");

    // Interval subdivision
    Stoer_Mult(rhs, a, b);
  }

  class CurveTroesch : public Curve
  {
  public:
    using Curve::Curve;

    VectorD2 operator()(FP_Type t)
    {
      VectorD2 y(2);
      y[0] = t;
      y[1] = 1;

      return y;
    }
  };

  std::pair<std::vector<FP_Type>, VectorD2>
  Troesch_MS(TimeFunctor &rhs, size_t dim, FP_Type a, FP_Type b)
  {
    // Approximate to BVP solution
    CurveTroesch* eta = new CurveTroesch(2);
    assert((*eta)(0)[0] == 0);
    assert((*eta)(1)[0] == 1);

    // Interval subdivision
    std::vector<FP_Type> t = trajectory(a, b, rhs, eta, 2, false, 1e-3);

    if (!(t.front() == a && t.back() == b))
      throw std::domain_error("subdivision does not match interval boundaries");

    size_t m = t.size();
    std::cout << "Amount of intervals: " << m-1 << std::endl;

    // Plot trajectory
    std::ofstream output_file;
    GnuPlot Dat1("Troesch_trajectory.dat", output_file);

    for (auto &c : t)
      output_file << c << "\t" << (*eta)(c);
    Dat1.plot_with_lines(2, "linespoints");

    // Starting values for Newton iteration
    VectorD2 s0(m*dim);

    for (size_t i = 0; i < m; i++)
      {
        VectorD2 s_i = (*eta)(t.at(i));

        for (size_t k = 0; k < dim; k++)
          s0[k + i*dim] = s_i[k];
      }
    return std::make_pair(t, s0);
  }

  void Troesch()
  {
    // Troesch BVP, see Num. Math. 2, 7.4.3.8
    FP_Type a = 0.0;
    FP_Type b = 1.0;

    // Linear separated BVP
    MatrixD2 A = init_matrix(2, 2, {1, 0, 0, 0});
    MatrixD2 B = init_matrix(2, 2, {0, 0, 1, 0});
    VectorD2 c = init_vector(2, {0, 1});
    BC_Linear r(A, B, c);

    // Left solution
    VectorD2 s(2);
    s[0] = 0;
    s[1] = 0.05;

    // 1. single shooting, external differentiation
    std::cout << "Single shooting (Troesch, ext. diff.)"
              << std::endl;

    std_tWrapper rhs(RHS_Troesch<VectorD2>, 2);
    SF_External<DOPRI54> M(rhs, true, 1e-3, 1e-4);

    SingleShooting F(M, a, b, r);
    Newton N(F, 2);
    N.iterate(s);

    // 2. single shooting, automatic differentation
    std::cout << "Single shooting (Troesch, aut. diff.)"
              << std::endl;

    FAD_tWrapper rhs_ad(RHS_Troesch<VectorAD>, 2);
    SF_Automatic<DOPRI54> M_ad(rhs_ad, true, 1e-3, 1e-4);

    SingleShooting F_ad(M_ad, a, b, r);
    Newton N_ad(F_ad, 2);
    N_ad.iterate(s);

    // 3. multiple shooting, external differentation
    auto T = Troesch_MS(rhs, 2, a, b);
    size_t m = T.first.size();
    std::cout << "Multiple shooting (Troesch, ext. diff.)"
              << std::endl;

    MultipleShooting G(M, T.first, r);
    Newton O(G, 2*m);
    O.iterate(T.second);

    // 4. multiple shooting, automatic differentiation
    std::cout << "Multiple shooting (Troesch, aut. diff.)"
              << std::endl;

    MultipleShooting G_ad(M_ad, T.first, r);
    Newton O_ad(G_ad, 2*m);
    O_ad.iterate(T.second);
  }
}

#endif // BVP_EXTERNAL_H
