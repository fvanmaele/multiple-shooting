#ifndef FORWARD_AD_H
#define FORWARD_AD_H
// Based on examples from:
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/ad_example.cpp
//
// For general information, see:
// https://software.sandia.gov/SESS/past_seminars/111307_Phipps.html

#include <cassert>
#include <vector>

#include <deal.II/lac/full_matrix.h>
#include <Sacado.hpp>

#include "functor.h"
#include "types.h"

// Forward AD class (dynamic allocation)
typedef Sacado::Fad::DFad<FP_Type> FAD_Number;

struct FAD_Functor
{
  virtual std::vector<FAD_Number>
  operator()(const std::vector<FAD_Number> &x) = 0;

  virtual ~FAD_Functor()
  {}
};

// Class representing a function f: R^m -> R^n, with support
// for evaluation and automatic (forward) differentation.
class FAD_Wrapper
{
public:
  FAD_Wrapper(FAD_Functor &_f, size_t _m, size_t _n) :
    f(_f), FAD_x(_m), FAD_y(_n), m(_m), n(_n)
  {}

  // Instantiate AD template classes and functions
  void init(const dealii::Vector<FP_Type> &x)
  {
    // Ensure argument length matches arguments
    assert(x.size() == m);

    for (size_t i = 0; i < x.size(); i++)
      {
        FAD_x.at(i) = x[i];
        FAD_x.at(i).diff(i, m);
      }

    FAD_y = f(FAD_x);
    assert(FAD_y.size() == n);

    FAD_initialized = true;
  }

  // Evaluate function
  dealii::Vector<FP_Type> value() const
  {
    if (!FAD_initialized)
      {
        throw std::invalid_argument("FAD must be initialized");
      }
    dealii::Vector<FP_Type> result(n);

    for (size_t i = 0; i < n; i++)
      {
        result(i) = FAD_y.at(i).val();
      }
    return result;
  }

  // Evaluate Jacobian
  dealii::FullMatrix<FP_Type> jacobian()
  {
    if (!FAD_initialized)
      {
        throw std::invalid_argument("FAD must be initialized");
      }
    dealii::FullMatrix<FP_Type> J(m, n);

    for (size_t i = 0; i < n; i++)
      {
        if (FAD_y.at(i).hasFastAccess())
            for (size_t j = 0; j < m; j++)
              {
                J.set(i, j, FAD_y.at(i).fastAccessDx(j));
              }
        else
            for (size_t j = 0; j < m; j++)
              {
                J.set(i, j, FAD_y.at(i).dx(j));
              }
      }
    return J;
  }

private:
  // Function definition
  FAD_Functor &f;

  // Save FAD input values in std::vector.
  // See: https://github.com/dealii/dealii/issues/6940
  std::vector<FAD_Number> FAD_x;
  std::vector<FAD_Number> FAD_y;

  // Markers
  bool FAD_initialized;
  const size_t m;
  const size_t n;
};

// This class assumes that the comparison of two vectors
// is less expensive than the re-evaluation of functions
// or their Jacobians. (cf. Sacado "FastAccess"?)
class Function_AD : public DivFunctor
{
public:
  Function_AD(FAD_Functor &f, size_t m, size_t n) :
    F(f, m, n), x_last(0)
  {}

  virtual dealii::Vector<FP_Type>
  value(const dealii::Vector<FP_Type> &x) override
  {
    if (!(x.size() == x_last.size() && x == x_last))
      {
        F.init(x);
        x_last = x;
      }
    return F.value();
  }

  virtual dealii::FullMatrix<FP_Type>
  jacobian(const dealii::Vector<FP_Type> &x) override
  {
    if (!(x.size() == x_last.size() && x == x_last))
      {
        F.init(x);
        x_last = x;
      }
    return F.jacobian();
  }

private:
  FAD_Wrapper F;
  dealii::Vector<FP_Type> x_last;
};

#endif // FORWARD_AD_H
