#ifndef VECTOR_OPERATORS_H
#define VECTOR_OPERATORS_H
#include <deal.II/lac/vector.h>
#include <algorithm>
#include "../base/number_type.h"

// The following operators are added for better readability
// of vector operations, at the cost of efficiency (additional copying
// of vectors is required).
dealii::Vector<NumberType> operator*(NumberType a, dealii::Vector<NumberType> rhs)
{
  rhs *= a;
  return rhs;
}

dealii::Vector<NumberType> operator*(dealii::Vector<NumberType> lhs, NumberType a)
{
  lhs *= a;
  return lhs;
}

dealii::Vector<NumberType> operator+(dealii::Vector<NumberType> lhs,
                                 const dealii::Vector<NumberType> &rhs)
{
  lhs += rhs;
  return lhs;
}

dealii::Vector<NumberType> operator-(dealii::Vector<NumberType> lhs,
                                 const dealii::Vector<NumberType> &rhs)
{
  lhs -= rhs;
  return lhs;
}

dealii::Vector<NumberType> std_to_Vector(const std::vector<NumberType> &v) {
  const size_t n = v.size();
  // create a dealii vector of size n and values set to 0
  dealii::Vector<NumberType> result(n);

  // create index vector 0...n-1
  std::vector<dealii::types::global_dof_index> indices(n, 0);
  std::iota(indices.begin(), indices.end(), 0);

  // add values stored in v to corresponding values in result
  result.add<NumberType>(indices, v);
  return result;
}

#endif // VECTOR_OPERATORS_H
