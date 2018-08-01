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

// C++17: std::optional for time parameter t
struct FAD_TdVecField
{
  virtual std::vector<FAD_Number>
  operator()(FAD_Number t, const std::vector<FAD_Number> &x) = 0;

  virtual ~FAD_TdVecField() = default;
};

// Class representing the right-hand side
//    f: I x R^d -> R^d
//
// of an IVP, supporting evaluation and automatic differentation.
template <size_t dim>
class FAD_Wrapper
{
public:
  FAD_Wrapper(FAD_TdVecField &_f) :
    f(_f), FAD_t(), FAD_u(dim), FAD_y(dim)
  {}

  // Instantiate AD template classes and functions
  void init(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == dim);

    // Time parameter (passive variable)
    FAD_t = t;

    // Analytic derivative with respect to u
    for (size_t i = 0; i < u.size(); i++)
      {
        FAD_u.at(i) = u[i];
        FAD_u.at(i).diff(i, dim);
      }

    FAD_y = f(FAD_t, FAD_u);
    FAD_initialized = true;
  }

  // Evaluate function
  dealii::Vector<FP_Type> value() const
  {
    assert(FAD_y.size() == dim);

    if (!FAD_initialized)
      {
        throw std::invalid_argument("FAD must be initialized");
      }
    dealii::Vector<FP_Type> result(dim);

    for (size_t i = 0; i < dim; i++)
      {
        result(i) = FAD_y.at(i).val();
      }
    return result;
  }

  // Evaluate partial derivatives with respect to u
  dealii::FullMatrix<FP_Type> nabla_u() const
  {
    assert(FAD_y.size() == dim);

    if (!FAD_initialized)
      {
        throw std::invalid_argument("FAD must be initialized");
      }
    dealii::FullMatrix<FP_Type> J(dim, dim);

    for (size_t i = 0; i < dim; i++)
      {
        if (FAD_y.at(i).hasFastAccess())
            for (size_t j = 0; j < dim; j++)
              {
                J.set(i, j, FAD_y.at(i).fastAccessDx(j));
              }
        else
            for (size_t j = 0; j < dim; j++)
              {
                J.set(i, j, FAD_y.at(i).dx(j));
              }
      }
    return J;
  }

private:
  // Function definition
  FAD_TdVecField &f;
  FAD_Number FAD_t;

  // Save FAD input values in std::vector.
  // https://github.com/dealii/dealii/issues/6940
  std::vector<FAD_Number> FAD_u;
  std::vector<FAD_Number> FAD_y;

  // Markers
  bool FAD_initialized;
};

template <size_t dim>
class TimeFunctor_AD : public TimeDivFunctor
{
public:
  TimeFunctor_AD(FAD_TdVecField &f) :
    F(f) // once differentiable
  {}

  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    F.init(t, u);
    return F.value();
  }

  virtual dealii::FullMatrix<FP_Type>
  nabla_u(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    F.init(t, u);
    return F.nabla_u();
  }

private:
  FAD_Wrapper<dim> F;
};

#endif // FORWARD_AD_H
