#ifndef FUNCTOR_H
#define FUNCTOR_H
// Based on examples from:
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/ad_example.cpp
//
// For general information, see:
// https://software.sandia.gov/SESS/past_seminars/111307_Phipps.html

#include <cassert>
#include <vector>

#include <Sacado.hpp>

#include "../base/functor.h"
#include "../base/types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

// Forward AD class (dynamic allocation)
typedef Sacado::Fad::DFad<FP_Type> FAD_Number;

struct FAD_Functor
{
  virtual std::vector<FAD_Number>
  operator()(const std::vector<FAD_Number> &x) = 0;
};

// Class representing a scalar function f: R^m -> R, with support
// for evaluation and automatic (forward) differentation.
template <size_t m, size_t n>
class FAD_Wrapper
{
public:
  FAD_Wrapper(FAD_Functor &_f) :
    f(_f), FAD_x(m), FAD_y(n)
  {}

  // Instantiate AD template classes and functions
  void init(const dealii::Vector<FP_Type> &x)
  {
    // Ensure argument length matches template
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
  dealii::LAPACKFullMatrix<FP_Type> jacobian()
  {
    if (!FAD_initialized)
      {
        throw std::invalid_argument("FAD must be initialized");
      }
    dealii::LAPACKFullMatrix<FP_Type> J(m, n);

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
};


#endif // FUNCTOR_H
