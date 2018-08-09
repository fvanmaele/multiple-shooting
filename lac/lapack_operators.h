#ifndef LAPACK_OPERATORS_H
#define LAPACK_OPERATORS_H

#include <iostream>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>

template class dealii::LAPACKFullMatrix<long double>;

#include "../base/types.h"

// LAPACKFullMatrix operators

dealii::Vector<FP_Type>
operator*(const dealii::LAPACKFullMatrix<FP_Type> &A,
          const dealii::Vector<FP_Type> &v)
{
  assert(A.n() == v.size());
  dealii::Vector<FP_Type> result(v.size());
  A.vmult(result, v);

  return result;
}

dealii::LAPACKFullMatrix<FP_Type>
operator*(const dealii::LAPACKFullMatrix<FP_Type> &A,
          const dealii::LAPACKFullMatrix<FP_Type> &B)
{
  assert(A.n() == B.m());
  dealii::LAPACKFullMatrix<FP_Type> result(A.m(), B.n());
  A.mmult(result, B);

  return result;
}

dealii::LAPACKFullMatrix<FP_Type>
operator*(FP_Type a, dealii::LAPACKFullMatrix<FP_Type> rhs)
{
  rhs *= a;
  return rhs;
}

dealii::LAPACKFullMatrix<FP_Type>
operator*(dealii::LAPACKFullMatrix<FP_Type> lhs, FP_Type a)
{
  lhs *= a;
  return lhs;
}

dealii::LAPACKFullMatrix<FP_Type>
operator+(dealii::LAPACKFullMatrix<FP_Type> result,
          const dealii::LAPACKFullMatrix<FP_Type> &B)
{
  assert(result.n() == B.n());
  assert(result.m() == B.m());
  result.add(1, B);

  return result;
}

std::ostream&
operator<<(std::ostream &os, const dealii::LAPACKFullMatrix<FP_Type> &A)
{
  A.print_formatted(os, 3, true, 0, "0");
  return os;
}

#endif // LAPACK_OPERATORS_H
