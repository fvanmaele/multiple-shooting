#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H

#include <iostream>

#include "../base/types.h"
#include "lac_types.h"

// FullMatrix operators

dealii::Vector<FP_Type>
operator*(const dealii::FullMatrix<FP_Type> &A,
          const dealii::Vector<FP_Type> &v)
{
  assert(A.n() == v.size());
  dealii::Vector<FP_Type> result(v.size());
  A.vmult(result, v);

  return result;
}

dealii::FullMatrix<FP_Type>
operator*(const dealii::FullMatrix<FP_Type> &A,
          const dealii::FullMatrix<FP_Type> &B)
{
  assert(A.n() == B.m());
  dealii::FullMatrix<FP_Type> result(A.m(), B.n());
  A.mmult(result, B);

  return result;
}

dealii::FullMatrix<FP_Type>
operator*(FP_Type a, dealii::FullMatrix<FP_Type> rhs)
{
  rhs *= a;
  return rhs;
}

dealii::FullMatrix<FP_Type>
operator*(dealii::FullMatrix<FP_Type> lhs, FP_Type a)
{
  lhs *= a;
  return lhs;
}

dealii::FullMatrix<FP_Type>
operator+(dealii::FullMatrix<FP_Type> result,
          const dealii::FullMatrix<FP_Type> &B)
{
  assert(result.n() == B.n());
  assert(result.m() == B.m());
  result.add(1, B);

  return result;
}

std::ostream&
operator<<(std::ostream &os, const dealii::FullMatrix<FP_Type> &A)
{
  A.print_formatted(os, 3, true, 0, "0");
  return os;
}

dealii::FullMatrix<FP_Type>
init_matrix(size_t m, size_t n, std::initializer_list<FP_Type> list)
{
  // Check range
  assert(list.size() == m*n);

  dealii::FullMatrix<FP_Type> temp(m, n);
  size_t i = 0;
  size_t j = 0;

  // Fill values row-by-row
  for (auto it = list.begin(); it != list.end(); it++)
    {
      temp.set(i, j, *it);
      j++;

      if (j % m == 0)
        {
          i++;
          j=0;
        }
    }
  return temp;
}

#endif // MATRIX_OPERATORS_H
