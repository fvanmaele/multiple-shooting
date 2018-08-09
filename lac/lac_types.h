#ifndef LAC_TYPES_H
#define LAC_TYPES_H

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector.templates.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/full_matrix.templates.h>

#include "../base/types.h"

typedef dealii::Vector<FP_Type> VectorD2;
typedef dealii::FullMatrix<FP_Type> MatrixD2;

// Instantiations for extended precision types
// https://www.dealii.org/9.0.0/doxygen/deal.II/Instantiations.html
template class dealii::Vector<long double>;
template class dealii::FullMatrix<long double>;

class Curve
{
public:
  virtual VectorD2
  operator()(FP_Type t) = 0;

  virtual ~Curve() = default;
};

class TimeFunctor
{
public:
  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) = 0;

  virtual ~TimeFunctor() = default;
};

class TimeDivFunctor : public TimeFunctor
{
public:
  virtual MatrixD2
  diff(FP_Type t, const VectorD2 &u) = 0;

  virtual ~TimeDivFunctor() = default;
};

class Functor
{
public:
  virtual VectorD2
  operator()(const VectorD2 &x) = 0;

  virtual ~Functor() = default;
};

class DivFunctor : public Functor
{
public:
  virtual MatrixD2
  diff(const VectorD2 &x) = 0;

  virtual ~DivFunctor() = default;
};

template <typename Callable>
class std_cWrapper : public Functor
{
public:
  std_cWrapper(Callable _f, size_t _dim) :
    f(_f), dim(_dim)
  {}

  virtual dealii::Vector<FP_Type>
  operator()(const dealii::Vector<FP_Type> &u) override
  {
    assert(u.size() == dim);
    return f(u);
  }

private:
  Callable f;
  size_t dim;
};

template <typename Callable>
class std_tWrapper : public TimeFunctor
{
public:
  std_tWrapper(Callable _f, size_t _dim) :
    f(_f), dim(_dim)
  {}

  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    assert(u.size() == dim);
    return f(t, u);
  }

private:
  Callable f;
  size_t dim;
};

#endif // LAC_TYPES_H
