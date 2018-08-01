#ifndef FUNCTOR_H
#define FUNCTOR_H

#include "types.h"
#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"

// XXX: Expand functors to include dimension of domain/image
// (cf. dealii::Function)
class Functor
{
public:
  virtual dealii::Vector<FP_Type>
  value(const dealii::Vector<FP_Type> &s) = 0;

  virtual ~Functor() = default;
};

class DivFunctor : public Functor
{
public:
  virtual dealii::FullMatrix<FP_Type>
  jacobian(const dealii::Vector<FP_Type> &s) = 0;

  virtual ~DivFunctor() = default;
};

// Note that there is no distinction between scalar and vector
// functions; both are represented as a dealii vector (of size 1
// in the scalar case).
class TimeFunctor
{
public:
  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u) = 0;

  virtual ~TimeFunctor() = default;
};

class TimeDivFunctor : public TimeFunctor
{
public:
  // Note this method represents the Jacobian with respect to
  // the components of u.
  virtual dealii::FullMatrix<FP_Type>
  nabla_u(FP_Type t, const dealii::Vector<FP_Type> &u) = 0;

  virtual ~TimeDivFunctor() = default;
};

class MatrixFunctor
{
  virtual dealii::FullMatrix<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u) = 0;
};

// Column-by-column RHS of variational equation
class FundMatrixFunctor : public TimeFunctor
{
public:
  FundMatrixFunctor(TimeDivFunctor _f, dealii::Vector<FP_Type> _phi) :
    f(_f), phi(_phi)
  {}

  virtual dealii::Vector<FP_Type>
  value(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    // matrix-vector product (matrix_operators.h)
    return f.nabla_u(t, u) * phi;
  }

private:
  TimeDivFunctor &f;
  dealii::Vector<FP_Type> phi;
};

#endif // FUNCTOR_H
