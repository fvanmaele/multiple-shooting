#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include <cmath>
#include <fstream>
#include <limits>

#include "../base/types.h"
#include "../bvp/linear.h"

// RHS for the initial value problem:
//    w'' = 1.5 * w^2
// resp. the system
//    (u1, u2)' = (u2, 1.5 * u1^2)
class RHS_Stoer : public TimeFunctor
{
public:
  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == 2);
    dealii::Vector<FP_Type> result(2);

    result[0] = u[1];
    result[1] = 1.5 * std::pow(u[0], 2);

    return result;
  }
};

class RHS_Troesch : public TimeFunctor
{
public:
  RHS_Troesch(FP_Type _Lambda) :
    Lambda(_Lambda)
  {}

  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    dealii::Vector<FP_Type> result(2);
    result[0] = u[1];
    result[1] = Lambda * std::sinh(Lambda * u[0]);

    return result;
  }

private:
  FP_Type Lambda;
};

// Stoer, Bulirsch, Num. Math 2, pp.192 (problem of 2nd order)
void Test_Stoer()
{
  RHS_Stoer f;
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

  // Single-shooting method for separated BVP
  SimpleBVP BVP(f, a, b, c);
  BVP.single_shooting(start);

  // Create plot
  std::ofstream output_file;
  output_file.open("bvp_sval.dat");
  assert(output_file.is_open());

  std::vector<dealii::Vector<FP_Type> > range;
  s[0] = 4;
  s[1] = -100;

  for (size_t i = 0; i < 100; i++)
    {
      s[1] += 1;
      range.emplace_back(s);
    }
  BVP.shooting_graph(2, range, output_file);

  std::system("gnuplot -p -e \"plot 'bvp_sval.dat' using 2:4 with lines\"");
}

// Randwertaufgabe von Troesch, siehe Num. Math. 2, 7.4.3.8
void Test_Troesch()
{
  RHS_Troesch f(5);
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

  // Single-shooting method for separated BVP
  SimpleBVP BVP(f, a, b, c);
  BVP.single_shooting(start);
}

#endif // BVP_EXTERNAL_H