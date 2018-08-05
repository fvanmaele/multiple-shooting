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
  Newton(Callable _f, size_t _dim, FP_Type _TOL = 1e-8,
         bool _ssc = true, size_t _ssc_lim = 20)
    :
      f(_f), dim(_dim), y(_dim), y_norm(-1), TOL(_TOL),
      steps(0), ssc_steps(0), ssc_lim(_ssc_lim), ssc(_ssc)
  {}

  dealii::Vector<FP_Type>
  step(const dealii::FullMatrix<FP_Type> &J,
       const dealii::Vector<FP_Type> &x)
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

    if (ssc)
      {
        // Candidate for next step x_{k+1}
        dealii::Vector<FP_Type> x_next(x.size());

        // j ist not guaranteed to be bounded (see Remark 4.2.4)
        for (size_t j = 0; j < ssc_lim; j++)
          {
            // Compute new value of x
            x_next = x - std::pow(0.5, j)*d;

            if (f(x_next).l2_norm() < y_norm)
              break;
            else
              ssc_steps++;
          }
        return x_next;
      }
    else
      return x - d;
  }

  dealii::Vector<FP_Type>
  iterate(dealii::Vector<FP_Type> s, size_t step_limit = 50)
  {
    static_assert(std::is_base_of<DivFunctor, Callable>::value,
                  "function is not differentiable");
    steps = 0;
    ssc_steps = 0;

    for (size_t k = 0; k < step_limit; k++)
      {
        // Perform step of (quasi-)Newton method
        s = step(f.diff(s), s);
        steps++;

        if (y_norm < TOL)
          {
            if (y_norm < 0)
              throw std::invalid_argument("norm must be positive");
            if (ssc)
              std::cout << "Solution of F: (" << steps << " steps, "
                        << ssc_steps << " ssc) " << s;
            else
              std::cout << "Solution of F: (" << steps << " steps) " << s;
            return s;
          }
      }
    std::cerr << "Warning: step limit exceeded" << std::endl;
    return s;
  }

private:
  Callable f;
  size_t dim;
  dealii::Vector<FP_Type> y;
  FP_Type y_norm;
  FP_Type TOL;

  // Step-size control
  size_t steps, ssc_steps, ssc_lim;
  bool ssc;
};

#endif // NEWTON_H
