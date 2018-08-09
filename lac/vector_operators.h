#ifndef VECTOR_OPERATORS_H
#define VECTOR_OPERATORS_H
#include <algorithm>
#include <cmath>
#include <iostream>

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

std::ostream&
operator<<(std::ostream &os, const std::vector<FP_Type> &v)
{
  for (size_t i = 0; i < v.size(); i++)
    {
      if (i == v.size()-1)
        os << v.at(i);
      else
        os << v.at(i) << "\t";
    }
  os << std::endl;
  return os;
}

dealii::Vector<FP_Type>
init_vector(size_t dim, std::initializer_list<FP_Type> list)
{
  // Check range
  assert(dim == list.size());

  dealii::Vector<FP_Type> temp(dim);
  size_t i = 0;

  for (auto it = list.begin(); it != list.end(); it++)
    {
      temp(i) = *it;
      i++;
    }
  return temp;
}

#endif // VECTOR_OPERATORS_H
