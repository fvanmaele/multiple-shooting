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

// Class representing the function
//    f: I x R^d -> R^d
//
// supporting evaluation and automatic differentation.
template <typename Callable>
class FAD_tWrapper : public TimeDivFunctor
{
public:
  FAD_tWrapper(Callable _f, size_t _dim) :
    f(_f), dim(_dim), f_init(false), FAD_t(0), FAD_u(_dim), FAD_y(_dim)
  {}

  void init(FP_Type t, const VectorD2 &u)
  {
    assert(u.size() == dim);

    // Time parameter (passive variable)
    FAD_t = t;

    // Analytic derivative with respect to u
    for (size_t i = 0; i < dim; i++)
      {
        FAD_u.at(i) = u[i];
        FAD_u.at(i).diff(i, dim);
      }

    FAD_y = f(FAD_t, FAD_u);
    f_init = true;
  }

  // Evaluate function
  VectorD2 value() const
  {
    if (!f_init)
      throw std::invalid_argument("FAD must be initialized");
    VectorD2 y(dim);

    for (size_t i = 0; i < FAD_y.size(); i++)
      y(i) = FAD_y.at(i).val();

    return y;
  }

  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) override
  {
    init(t, u);
    return value();
  }

  // Evaluate partial derivatives with respect to u
  MatrixD2 diff() const
  {
    if (!f_init)
      throw std::invalid_argument("FAD must be initialized");
    MatrixD2 J(dim, dim);

    for (size_t i = 0; i < dim; i++)
      if (FAD_y.at(i).hasFastAccess())
        for (size_t j = 0; j < dim; j++)
          J.set(i, j, FAD_y.at(i).fastAccessDx(j));
      else
        for (size_t j = 0; j < dim; j++)
          J.set(i, j, FAD_y.at(i).dx(j));

    return J;
  }

  virtual MatrixD2
  diff(FP_Type t, const VectorD2 &u) override
  {
    init(t, u);
    return diff();
  }

private:
  Callable f;
  size_t dim;
  bool f_init;

  NumberAD FAD_t;
  VectorAD FAD_u;
  VectorAD FAD_y;
};

// Class which fixes t in f(t, u).
template <typename Callable>
class FAD_cWrapper : public DivFunctor
{
public:
  FAD_cWrapper(Callable f, size_t dim, FP_Type _t = 0) :
    F(f, dim), t(_t)
  {}

  virtual VectorD2
  operator()(const VectorD2 &u) override
  {
    F.init(t, u);
    return F.value();
  }

  virtual MatrixD2
  diff(const VectorD2 &u) override
  {
    F.init(t, u);
    return F.diff();
  }

private:
  FAD_tWrapper<Callable> F;
  FP_Type t;
};

#endif // FORWARD_AD_H
