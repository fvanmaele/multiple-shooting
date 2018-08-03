#ifndef NEWTON_H
#define NEWTON_H

#include <deal.II/lac/lapack_full_matrix.h>

#include "../base/types.h"
#include "../base/forward_ad.h"
#include "../lac/lac_types.h"
#include "../lac/vector_operators.h"

// Solve nonlinear root finding problem:
//    f(x) = 0,   f: R^d -> R^d
template <typename Callable>
class Newton
{
public:
  // Step-size control requires that f is callable, since
  // new values of y are computed for varying x.
  Newton(Callable _f, size_t _dim) :
    f(_f), y(_dim), y_norm(-1), dim(_dim)
  {}

  dealii::Vector<FP_Type>
  step(const dealii::FullMatrix<FP_Type> &J,
       const dealii::Vector<FP_Type> &x,
       bool step_size_control = true,
       size_t ssc_limit = 20)
  {
    assert(x.size() == dim);
    y = f(x);
    y_norm = y.l2_norm();

    // The Jacobian J is taken argument to allow using
    // this function for both Newton and quasi-Newton methods.
    dealii::LAPACKFullMatrix<FP_Type> Jacobian(y.size(), y.size());
    Jacobian = J;

    // Solve the linear system. This is done in-place;
    // we save a copy of y for step-size control.
    Jacobian.compute_lu_factorization();
    dealii::Vector<FP_Type> d(y);
    Jacobian.solve(d);
    if (step_size_control)
      {
        // Candidate for next step x_{k+1}
        dealii::Vector<FP_Type> x_next(x.size());

        // j ist not guaranteed to be bounded (see Remark 4.2.4)
        for (size_t j = 0; j < ssc_limit; j++)
          {
            // Compute new value of x
            x_next = x - std::pow(0.5, j)*d;

            if (f(x_next).l2_norm() < y_norm)
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

  dealii::Vector<FP_Type>
  iterate(dealii::Vector<FP_Type> s, size_t step_limit = 50)
  {
    static_assert(std::is_base_of<DivFunctor, Callable>::value,
                  "function is not differentiable");
    size_t steps = 0;

    for (size_t k = 0; k < step_limit; k++)
      {
        // Perform step of (quasi-)Newton method
        s = step(f.diff(s), s, true);
        steps++;

        if (stopping_criterion(1e-8))
          {
            std::cout << "Solution of F: (" << steps << " steps) " << s;
            return s;
          }
      }
    std::cerr << "Warning: step limit exceeded" << std::endl;
    return s;
  }

private:
  Callable f;
  dealii::Vector<FP_Type> y;
  FP_Type y_norm;
  size_t dim;
};

#endif // NEWTON_H
