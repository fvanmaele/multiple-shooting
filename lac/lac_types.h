#ifndef LAC_TYPES_H
#define LAC_TYPES_H

#include <functional>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include "../base/types.h"

typedef dealii::Vector<FP_Type> VectorD2;
typedef std::function<VectorD2(FP_Type, VectorD2)> tVecField;
typedef std::function<VectorD2(VectorD2)> cVecField;

typedef dealii::FullMatrix<FP_Type> MatrixD2;
typedef std::function<MatrixD2(FP_Type, VectorD2)> tMatField;
typedef std::function<MatrixD2(VectorD2)> cMatField;

class TimeDivFunctor
{
public:
  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) = 0;

  virtual MatrixD2
  diff(FP_Type t, const VectorD2 &u) = 0;

  virtual ~TimeDivFunctor() = default;
};

class DivFunctor
{
public:
  virtual VectorD2
  operator()(const VectorD2 &x) = 0;

  virtual MatrixD2
  diff(const VectorD2 &x) = 0;

  virtual ~DivFunctor() = default;
};

#endif // LAC_TYPES_H
