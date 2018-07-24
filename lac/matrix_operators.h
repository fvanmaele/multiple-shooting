#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "../base/types.h"

// LAPACK matrices are stored in transposed order, requiring
// a roundabout to FullMatrix.
template <size_t dim>
dealii::LAPACKFullMatrix<FP_Type>
MFA(size_t rows, size_t cols, const std::array<FP_Type, dim> &entries)
{
  assert(dim == rows * cols);
  dealii::LAPACKFullMatrix<FP_Type> A(rows, cols);

  A = dealii::FullMatrix<FP_Type>(rows, cols, entries.data());
  return A;
}

// Matrix/vector multiplication
dealii::Vector<FP_Type>
operator*(dealii::LAPACKFullMatrix<FP_Type> &A,
          const dealii::Vector<FP_Type> &x)
{
  assert(false); // to implement
}

#endif // MATRIX_OPERATORS_H
