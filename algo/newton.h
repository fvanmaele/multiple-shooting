#ifndef NEWTON_H
#define NEWTON_H

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>

template class dealii::LAPACKFullMatrix<long double>;

#include "../base/types.h"
#include "../base/forward_ad.h"
#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

// Solve nonlinear root finding problem:
//    f(x) = 0,   f: R^d -> R^d
template <typename Callable>
class Newton
{
public:
  // Step-size control requires that f is callable, since
  // new values of y are computed for varying x.
  Newton(Callable _f, size_t _dim, FP_Type _TOL = 1e-6,
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

    // The Jacobian J is taken as argument to allow using
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
  iterate(dealii::Vector<FP_Type> x, size_t step_limit = 25)
  {
    static_assert(std::is_base_of<DivFunctor, Callable>::value,
                  "function is not differentiable");
    steps = 0;
    ssc_steps = 0;

    for (size_t k = 0; k < step_limit; k++)
      {
        // Perform step of (quasi-)Newton method
        x = step(f.diff(x), x);
        steps++;

        if (y_norm < TOL)
          {
            if (y_norm < 0)
              throw std::invalid_argument("norm must be positive");
            if (ssc)
              std::cout << "Solution of F: (" << steps << " steps, "
                        << ssc_steps << " ssc) " << x;
            else
              std::cout << "Solution of F: (" << steps << " steps) " << x;
            return x;
          }
      }
    std::cerr << "Warning: step limit exceeded" << std::endl;
    return x;
  }

  dealii::Vector<FP_Type>
  iterate_broyden(dealii::Vector<FP_Type> x, size_t step_limit = 50)
  {
    std::cerr << "Using Broyden method" << std::endl;
    static_assert(std::is_base_of<DivFunctor, Callable>::value,
                  "function is not differentiable");
    steps = 0;
    ssc_steps = 0;

    dealii::FullMatrix<FP_Type> J = f.diff(x);
    dealii::FullMatrix<FP_Type> J_prev(J);
    dealii::Vector<FP_Type> x_prev(x);
    x = step(J, x);

    for (size_t k = 1; k < step_limit; k++)
      {
        // Full computation of Jacobian every 4 steps
        if (k % 4 == 0)
          {
            J = f.diff(x);
            x = step(J, x);
          }
        else
          {
            dealii::Vector<FP_Type> p = x - x_prev;
            dealii::Vector<FP_Type> q = f(x) - f(x_prev);

            dealii::FullMatrix<FP_Type> V(p.size(), p.size());
            V.outer_product(q - J_prev * p, p);

            // For larger J, use Sherman-Morrison formula to update J^{-1} directly
            J = J_prev + 1. / p.norm_sqr() * V;
            x = step(J, x);
          }
        steps++;

        if (y_norm < TOL)
          {
            if (y_norm < 0)
              throw std::invalid_argument("norm must be positive");
            if (ssc)
              std::cout << "Solution of F: (" << steps << " steps, "
                        << ssc_steps << " ssc) " << x;
            else
              std::cout << "Solution of F: (" << steps << " steps) " << x;
            return x;
          }
      }
    std::cerr << "Warning: step limit exceeded" << std::endl;
    return x;
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
