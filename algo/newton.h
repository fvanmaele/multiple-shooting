#ifndef NEWTON_H
#define NEWTON_H

#include <deal.II/lac/lapack_full_matrix.h>

#include "../base/functor.h"
#include "../base/types.h"
#include "../base/forward_ad.h"

// Solve nonlinear root finding problem:
//    f(x) = 0,   f: R^d -> R^d
class Newton
{
public:
  Newton(Functor &_F, size_t dim) :
    F(_F), y(dim), y_norm(-1)
  {}

  dealii::Vector<FP_Type>
  step(const dealii::FullMatrix<FP_Type> &Jacobian,
       const dealii::Vector<FP_Type> &x,
       bool step_size_control = true,
       size_t ssc_limit = 20)
  {
    y = F.value(x);
    y_norm = y.l2_norm();

    // The Jacobian J is taken argument to allow using
    // this function for both Newton and quasi-Newton methods.
    dealii::LAPACKFullMatrix<FP_Type> J(y.size(), y.size());
    J = Jacobian;

    // Solve the linear system. This is done in-place;
    // we save a copy of y for step-size control.
    J.compute_lu_factorization();
    dealii::Vector<FP_Type> d(y);
    J.solve(d);

    // This section requires that F is callable, since
    // new values of y are computed for varying x.
    if (step_size_control)
      {
        // Candidate for next step x_{k+1}
        dealii::Vector<FP_Type> x_next(x.size());

        // j ist not guaranteed to be bounded (see Remark 4.2.4)
        for (size_t j = 0; j < ssc_limit; j++)
          {
            // Compute new value of x
            x_next = x - std::pow(0.5, j)*d;

            if (F.value(x_next).l2_norm() < y_norm)
              break;
          }
        return x_next;
      }
    else
      return x - d;
  }

  bool stopping_criterion(FP_Type TOL) const
  {
    if (y_norm < 0)
      throw std::invalid_argument("no Newton step performed");

    return y_norm < TOL;
  }

private:
  Functor &F;
  dealii::Vector<FP_Type> y;
  FP_Type y_norm;
};

#endif // NEWTON_H
