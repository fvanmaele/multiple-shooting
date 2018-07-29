#ifndef LINEAR_H
#define LINEAR_H

#include <deal.II/lac/full_matrix.h>

#include "../base/functor.h"
#include "../base/types.h"

// Class to solve a linear BVP
//   y' = f(x, y),  A*y(a) + B*y(b) = c
class BVP
{
public:

private:
  RHS &f;
  dealii::FullMatrix<FP_Type> A;
  dealii::FullMatrix<FP_Type> B;
  FP_Type a, b, c;
};

#endif // LINEAR_H
