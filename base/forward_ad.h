#ifndef FORWARD_AD_H
#define FORWARD_AD_H
// Based on examples from:
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/ad_example.cpp
//
// For general information, see:
// https://software.sandia.gov/SESS/past_seminars/111307_Phipps.html

#include <cassert>
#include <vector>

#include <Sacado.hpp>

#include "types.h"
#include "../lac/lac_types.h"

// Forward AD class (dynamic allocation)
typedef Sacado::Fad::DFad<FP_Type> NumberAD;

// Save FAD input values in std::vector.
// https://github.com/dealii/dealii/issues/6940
typedef std::vector<NumberAD> VectorAD;

// Helper functions
void fad_differentiate(const std::vector<NumberAD> &FAD_y,
                       dealii::FullMatrix<FP_Type> &J)
{
  size_t dim = FAD_y.size();

  for (size_t i = 0; i < dim; i++)
    if (FAD_y.at(i).hasFastAccess())
      for (size_t j = 0; j < dim; j++)
        J.set(i, j, FAD_y.at(i).fastAccessDx(j));
    else
      for (size_t j = 0; j < dim; j++)
        J.set(i, j, FAD_y.at(i).dx(j));
}

void fad_evaluate(const std::vector<NumberAD> &FAD_y,
                  dealii::Vector<FP_Type> &y)
{
  for (size_t i = 0; i < FAD_y.size(); i++)
    y(i) = FAD_y.at(i).val();
}

void fad_set_vars(std::vector<NumberAD> &FAD_u,
                  const dealii::Vector<FP_Type> &u)
{
  size_t dim = FAD_u.size();

  for (size_t i = 0; i < dim; i++)
    {
      FAD_u.at(i) = u[i];
      FAD_u.at(i).diff(i, dim);
    }
}

// Class representing the function
//    f: R^d -> R^d
//
// supporting evaluation and automatic differentation.
template <typename Callable>
class FAD_cWrapper : public DivFunctor
{
public:
  FAD_cWrapper(Callable _f, size_t _dim) :
    f(_f), dim(_dim), FAD_u(_dim), FAD_y(_dim)
  {}

  // Instantiate AD template classes and functions
  void init(const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == dim);

    // Analytic derivative with respect to u
    fad_set_vars(FAD_u, u);

    FAD_y = f(FAD_u);
    FAD_initialized = true;
  }

  // Evaluate function
  dealii::Vector<FP_Type>
  value() const
  {
    if (!FAD_initialized)
      throw std::invalid_argument("FAD must be initialized");
    dealii::Vector<FP_Type> y(dim);

    fad_evaluate(FAD_y, y);
    return y;
  }

  virtual dealii::Vector<FP_Type>
  operator()(const dealii::Vector<FP_Type> &u) override
  {
    init(u);
    return value();
  }

  // Evaluate partial derivatives with respect to u
  dealii::FullMatrix<FP_Type>
  diff() const
  {
    if (!FAD_initialized)
      throw std::invalid_argument("FAD must be initialized");
    dealii::FullMatrix<FP_Type> J(dim, dim);

    fad_differentiate(FAD_y, J);
    return J;
  }

  virtual dealii::FullMatrix<FP_Type>
  diff(const dealii::Vector<FP_Type> &u) override
  {
    init(u);
    return diff();
  }

private:
  Callable f;
  size_t dim;
  std::vector<NumberAD> FAD_u;
  std::vector<NumberAD> FAD_y;
  bool FAD_initialized;
};

// Class representing the function
//    f: I x R^d -> R^d
//
// supporting evaluation and automatic differentation.
template <typename Callable>
class FAD_tWrapper : public TimeDivFunctor
{
public:
  FAD_tWrapper(Callable _f, size_t _dim) :
    f(_f), FAD_t(0), dim(_dim), FAD_u(_dim), FAD_y(_dim)
  {}

  void init(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == dim);

    // Time parameter (passive variable)
    FAD_t = t;

    // Analytic derivative with respect to u
    fad_set_vars(FAD_u, u);

    FAD_y = f(FAD_t, FAD_u);
    FAD_initialized = true;
  }

  // Evaluate function
  dealii::Vector<FP_Type>
  value() const
  {
    if (!FAD_initialized)
      throw std::invalid_argument("FAD must be initialized");
    dealii::Vector<FP_Type> y(dim);

    fad_evaluate(FAD_y, y);
    return y;
  }

  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    init(t, u);
    return value();
  }

  // Evaluate partial derivatives with respect to u
  dealii::FullMatrix<FP_Type>
  diff() const
  {
    if (!FAD_initialized)
      throw std::invalid_argument("FAD must be initialized");
    dealii::FullMatrix<FP_Type> J(dim, dim);

    fad_differentiate(FAD_y, J);
    return J;
  }

  virtual dealii::FullMatrix<FP_Type>
  diff(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    init(t, u);
    return diff();
  }

private:
  Callable f;
  NumberAD FAD_t;
  size_t dim;
  std::vector<NumberAD> FAD_u;
  std::vector<NumberAD> FAD_y;
  bool FAD_initialized;
};

#endif // FORWARD_AD_H
