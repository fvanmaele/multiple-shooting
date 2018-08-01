#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H

#include <deal.II/lac/full_matrix.h>

#include "../base/types.h"

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

template <size_t n, size_t m>
dealii::FullMatrix<FP_Type>
std_to_Matrix(const std::array<FP_Type, n*m> &v)
{
  return dealii::FullMatrix<FP_Type>(n, m, v.data());
}

#endif // MATRIX_OPERATORS_H
