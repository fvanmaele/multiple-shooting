#ifndef VECTOR_OPERATORS_H
#define VECTOR_OPERATORS_H
#include <algorithm>

#include "../base/types.h"
#include "lac_types.h"

// The following operators are added for better readability
// of vector operations, at the cost of efficiency (additional copying
// of vectors is required).
dealii::Vector<FP_Type>
operator*(const FP_Type a, dealii::Vector<FP_Type> rhs)
{
  rhs *= a;
  return rhs;
}

dealii::Vector<FP_Type>
operator*(dealii::Vector<FP_Type> lhs, const FP_Type a)
{
  lhs *= a;
  return lhs;
}

dealii::Vector<FP_Type>
operator/(dealii::Vector<FP_Type> lhs, const FP_Type a)
{
  lhs /= a;
  return lhs;
}

dealii::Vector<FP_Type>
operator+(dealii::Vector<FP_Type> lhs,
          const dealii::Vector<FP_Type> &rhs)
{
  lhs += rhs;
  return lhs;
}

dealii::Vector<FP_Type>
operator-(dealii::Vector<FP_Type> lhs,
          const dealii::Vector<FP_Type> &rhs)
{
  lhs -= rhs;
  return lhs;
}

#endif // VECTOR_OPERATORS_H
