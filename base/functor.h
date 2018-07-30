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
  value(const dealii::Vector<FP_Type> &x) = 0;

  virtual ~Functor()
  {}
};

class DivFunctor : public Functor
{
public:
  virtual dealii::FullMatrix<FP_Type>
  jacobian(const dealii::Vector<FP_Type> &s) = 0;

  virtual ~DivFunctor()
  {}
};

#endif // FUNCTOR_H
