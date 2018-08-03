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

dealii::Vector<FP_Type>
std_to_Vector(const std::vector<FP_Type> &v)
{
  const size_t n = v.size();
  // create a dealii vector of size n and values set to 0
  dealii::Vector<FP_Type> result(n);

  // create index vector 0...n-1
  std::vector<dealii::types::global_dof_index> indices(n, 0);
  std::iota(indices.begin(), indices.end(), 0);

  // add values stored in v to corresponding values in result
  result.add<FP_Type>(indices, v);
  return result;
}

#endif // VECTOR_OPERATORS_H
