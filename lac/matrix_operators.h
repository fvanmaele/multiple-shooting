#ifndef MATRIX_OPERATORS_H
#define MATRIX_OPERATORS_H
#include <deal.II/lac/full_matrix.h>
#include "../base/number_type.h"

NumberType matrix_element(const dealii::FullMatrix<NumberType> &A,
                      std::size_t i, std::size_t j)
{
  dealii::FullMatrix<NumberType>::Accessor element(&A, i, j);
  return element.value();
}

#endif // MATRIX_OPERATORS_H
