#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <vector>
#include <Sacado.hpp>
#include "../base/number_type.h"

// Based on examples from:
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/ad_example.cpp
// https://www.dealii.org/9.0.0/doxygen/deal.II/step_33.html#Trilinossolvers
//
// For general information, see:
// https://en.wikipedia.org/wiki/Automatic_differentiation
typedef Sacado::Fad::DFad<NumberType>  FAD_Number;
typedef Sacado::Rad::ADvar<NumberType> RAD_Number;

// Every component f_i, i = 0,..,n-1, of a function
//    f: Rm -> Rn   (m degrees of freedom)
//
// is stored as an FAD_Number object. The corresponding vector
// (including the values where the derivate is evaluated) must
// be defined separately. For example, the component:
//    f: R^2 -> R,  f(x, y) = 2x + cos(xy)
//
// is defined at (x, y) = (1, 2) by:
//    FAD_Number x = 1;
//    FAD_Number y = 2;
//    x.diff(0, 2); // 1st partial derivative (of 2)
//    y.diff(1, 2); // 2nd partial derivative
//    FAD_Number f = 2*x + cos(x*y);
//
// The return value is an m*n Matrix containing all evaluated
// partial derivatives (NumberType).
class FAD_Jacobian
{
public:
  // https://github.com/dealii/dealii/issues/6940
  template <size_t dim>
  dealii::FullMatrix<NumberType> diff(const dealii::Tensor<1, dim, FAD_Number> &f, size_t m)
  {
    size_t n = f.dimension;
    dealii::FullMatrix<NumberType> J(n, m);

    // Fill entries of Jacobi matrix
    for (size_t i = 0; i < n; i++)
      {
        for (size_t j = 0; j < m; j++)
          {
            J.set(i, j, f[i].fastAccessDx(j));
          }
      }
    return J;
  }
};

void Test_FAD()
{
  // (x,y) = (1,2)
  FAD_Number x = 1;
  FAD_Number y = 2;

  // f: R^2 -> R^3,
  x.diff(0, 2);
  y.diff(1, 2);

  dealii::Tensor<1, 3, FAD_Number> f;
  f[0] = pow(x, 3) * sin(y);
  f[1] = cos(x) * sin(y);
  f[2] = exp(x);

  // Output Jacobi matrix
  FAD_Jacobian Method;
  Method.diff<3>(f, 2).print_formatted(std::cout);
}

#endif // JACOBIAN_H
