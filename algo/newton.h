#ifndef NEWTON_H
#define NEWTON_H

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
      f(_f), dim(_dim), y_norm(-1), TOL(_TOL),
      steps(0), ssc_steps(0), ssc_lim(_ssc_lim), ssc(_ssc)
  {}

  // The Jacobian J is taken as argument to allow using
  // this function for both Newton and quasi-Newton methods.
  void step(MatrixD2 &J, VectorD2 &x)
  {
    assert(x.size() == dim);
    VectorD2 y = f(x);

    assert(y.size() == dim);
    y_norm = y.l2_norm();

    // Solve the linear system in-place.
    J.gauss_jordan();
    VectorD2 d = J * y;

    if (ssc)
      { // Candidate for next step x_{k+1}
        VectorD2 x_next(x.size());

        // j ist not guaranteed to be bounded (see Remark 4.2.4)
        for (size_t j = 0; j < ssc_lim; j++)
          {
            x_next = x - std::pow(0.5, j) * d;

            if (f(x_next).l2_norm() < y_norm)
              break;
            else
              ssc_steps++;
          }

        steps++;
        x = x_next;
      }
    else
      {
        steps++;
        x = x - d;
      }
  }

  // This method assumes f is differentiable, i.e. includes a diff() method.
  VectorD2 iterate(const VectorD2 &x0, size_t step_limit = 25)
  {
    steps = 0;
    ssc_steps = 0;
    VectorD2 x = x0;
    MatrixD2 J(dim, dim);

    for (size_t k = 0; k < step_limit; k++)
      {
        // Perform step of (quasi-)Newton method
        J = f.diff(x);
        step(J, x);

        if (y_norm < TOL)
          {
            if (y_norm < 0)
              throw std::invalid_argument("norm must be positive");

            if (ssc)
              {
                std::cout << "Solution of F: ("  << steps << " steps, "
                          << ssc_steps << " ssc) ";
                x.print(std::cout);
              }
            else
              {
                std::cout << "Solution of F: (" << steps << " steps) ";
                x.print(std::cout);
              }
            return x;
          }
      }

    std::cerr << "Warning: step limit exceeded" << std::endl;
    return x;
  }

  VectorD2 iterate_broyden(const VectorD2 &x0, size_t step_limit = 50)
  {
    VectorD2 x = x0;
    MatrixD2 J = f.diff(x);
    steps = 1;
    ssc_steps = 0;

    MatrixD2 J_prev(J);
    VectorD2 x_prev(x);
    step(J, x);

    for (size_t k = 1; k < step_limit; k++)
      {
        // Full computation of Jacobian every 4 steps
        if (k % 4 == 0)
          {
            J = f.diff(x);
            step(J, x);
          }
        else
          {
            // Rank-1 updates
            VectorD2 p = x - x_prev;
            VectorD2 q = f(x) - f(x_prev);
            MatrixD2 V(p.size(), p.size());

            // XXX: for large J, use Sherman-Morrison formula to update J^{-1} directly
            V.outer_product(q - J_prev * p, p);
            J = J_prev + 1. / p.norm_sqr() * V;
            step(J, x);
          }

        if (y_norm < TOL)
          {
            if (y_norm < 0)
              throw std::invalid_argument("norm must be positive");

            if (ssc)
              std::cout << "Solution of F (Broyden): (" << steps << " steps, "
                        << ssc_steps << " ssc) " << x;
            else
              std::cout << "Solution of F (Broyden): (" << steps << " steps) " << x;

            return x;
          }
      }

    std::cerr << "Warning: step limit exceeded" << std::endl;
    return x;
  }

private:
  Callable f;
  size_t dim;
  FP_Type y_norm;
  FP_Type TOL;

  // Step-size control
  size_t steps, ssc_steps, ssc_lim;
  bool ssc;
};

#endif // NEWTON_H
