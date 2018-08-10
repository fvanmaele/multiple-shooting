#ifndef LAC_TYPES_H
#define LAC_TYPES_H

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector.templates.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/full_matrix.templates.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_vector.templates.h>

#include "../base/types.h"

typedef dealii::Vector<FP_Type> VectorD2;
typedef dealii::BlockVector<FP_Type> BlockD2;
typedef dealii::FullMatrix<FP_Type> MatrixD2;

// Instantiations for extended precision types
// https://www.dealii.org/9.0.0/doxygen/deal.II/Instantiations.html
template class dealii::Vector<long double>;
template class dealii::FullMatrix<long double>;
template class dealii::BlockVector<long double>;

class Curve
{
public:
  Curve(size_t n) : dim(n)
  {}

  virtual ~Curve() = default;

  virtual VectorD2
  operator()(FP_Type t) = 0;

  size_t n_dim() const
  {
    return dim;
  }

private:
  size_t dim;
};

class Functor
{
public:
  Functor(size_t n) : dim(n)
  {}

  virtual ~Functor() = default;

  virtual VectorD2
  operator()(const VectorD2 &x) = 0;

  size_t n_dim() const
  {
    return dim;
  }

private:
  size_t dim;
};

class DivFunctor : public Functor
{
public:
  using Functor::Functor;

  virtual ~DivFunctor() = default;

  virtual MatrixD2
  diff(const VectorD2 &x) = 0;
};

class TimeFunctor
{
public:
  TimeFunctor(size_t n) : dim(n)
  {}

  virtual ~TimeFunctor() = default;

  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) = 0;

  size_t n_dim() const
  {
    return dim;
  }

private:
  size_t dim;
};

class TimeDivFunctor : public TimeFunctor
{
public:
  using TimeFunctor::TimeFunctor;

  virtual MatrixD2
  diff(FP_Type t, const VectorD2 &u) = 0;

  virtual ~TimeDivFunctor() = default;
};

template <typename Callable>
class std_cWrapper : public Functor
{
public:
  std_cWrapper(Callable _f, size_t dim) :
    Functor(dim), f(_f)
  {}

  virtual VectorD2
  operator()(const VectorD2 &u) override
  {
    return f(u);
  }

private:
  Callable f;
};

template <typename Callable>
class std_tWrapper : public TimeFunctor
{
public:
  std_tWrapper(Callable _f, size_t dim) :
    TimeFunctor(dim), f(_f)
  {}

  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) override
  {
    return f(t, u);
  }

private:
  Callable f;
};

#endif // LAC_TYPES_H
