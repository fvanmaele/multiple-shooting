#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "../base/types.h"

// LAPACK matrices are stored in transposed order, requiring
// a roundabout to FullMatrix.
dealii::LAPACKFullMatrix<FP_Type>
LAPACK_from_Array(size_t rows, size_t cols, const FP_Type *entries)
{
  dealii::LAPACKFullMatrix<FP_Type> A(rows, cols);

  A = dealii::FullMatrix<FP_Type>(rows, cols, entries);
  return A;
}

#endif // MATRIX_OPERATORS_H
