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

  template <typename DiffMethod>
  void Graph_Stoer(SimpleBVP<DiffMethod> &BVP, std::string filename)
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
  void Stoer()
  {
    std_tWrapper f(RHS_Stoer<VectorD2>, 2);
    FAD_tWrapper f_ad(RHS_Stoer<VectorAD>, 2);
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

    Graph_Stoer<SF_External>(BVP, "bvp_stoer_ed.dat");
    Graph_Stoer<SF_Automatic>(BVP_AD, "bvp_stoer_ad.dat");
  }

  class CurveLinear : public Curve
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

  // Troesch BVP, see Num. Math. 2, 7.4.3.8
  void Troesch()
  {
    std_tWrapper f(RHS_Troesch<VectorD2>, 2);
    FAD_tWrapper f_ad(RHS_Troesch<VectorAD>, 2);
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
    BVP_AD.single_shooting(start);

    // Interval subdivision for multiple shooting method
    CurveLinear* eta = new CurveLinear;
    assert((*eta)(0)[0] == 0);
    assert((*eta)(1)[0] == 1);

    std::vector<FP_Type> subint = trajectory(a, b, f, eta, 1e-2, 2, false);
    std::cout << "Amount of intervals: " << subint.size()-1 << std::endl;
    std::cout << subint;

    std::ofstream output_file;
    GnuPlot Dat1("Troesch_trajectory.dat", output_file);

    // Plot trajectory
    for (auto &c : subint)
      output_file << c << "\t" << (*eta)(c);
    Dat1.plot_with_lines(2, "linespoints");
  }
}

#endif // BVP_EXTERNAL_H
