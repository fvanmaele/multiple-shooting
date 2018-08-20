#ifndef NEWTON_H
#define NEWTON_H

#include "../base/types.h"
#include "../base/forward_ad.h"
#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

/*! \class Newton
 *  \brief Newton method with globalization and rank-1 updates.
 *
 * This class implements the Newton method for finding the root of a
 * non-linear equation \f$f:\mathbb{R}^d \rightarrow \mathbb{R}^d\f$.
 * Step size control is available as globalization strategy (Def 4.2.3)
 * in case a good intial guess of the root is not available.
 */
template <typename Callable>
class Newton
{
public:
  /*! \fn Newton
   *  \brief Constructor. Initialize the method with a function
   * \f$f:\mathbb{R}^d \rightarrow \mathbb{R}^d\f$ of dimension \f$d\f$.
   *
   * The TOL parameter specifies when a root \f$s\f$ is accepted. Step size control
   * is activated by default; as the value \f$j\f$ used within need not be
   * bounded, the maximum may be set here (Remark 4.2.4).
   */
  Newton(Callable _f, size_t _dim, FP_Type _TOL = 1e-6,
         bool _ssc = true, size_t _ssc_lim = 20)
    :
      f(_f), dim(_dim), y_norm(-1), TOL(_TOL),
      steps(0), ssc_steps(0), ssc_lim(_ssc_lim), ssc(_ssc)
  {}

  /*! \fn step
   *  \brief Perform a Newton step.
   *
   * The Jacobian \a inverse is specifically taken as argument, to allow use of
   * this function for both Newton and quasi-Newton methods.
   */
  void step(const MatrixD2 &J_inv, VectorD2 &x)
  {
    assert(x.size() == dim);
    VectorD2 y = f(x);

    assert(y.size() == dim);
    y_norm = y.l2_norm();

    // Solution of the linear system
    VectorD2 d = J_inv * y;

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

  /*! \fn iterate
   *  \brief Perform Newton steps until \f$\|f(s)\| < TOL\f$.
   *
   * As the Jacobian is computed in this function, we assume that \f$f\f$
   * is differentiable, i.e. has an available \c diff() method.
   *
   * The maximum amount of steps may be specified, defaulting to 25. In our context,
   * exceeding this limit has indicated either a program error or an unsuitably
   * chosen method. Should this occur, the function therefore exits with an exception.
   *
   * To solve the resulting linear systems, LU decomposition is used
   * via dealii and LAPACK. The Jacobian that results from the multiple
   * shooting method is sparse, but of small dimension in the problems we
   * consider. (In particular, the Thomas-Fermi problem with 20 subintervals
   * results in a 40 x 40 Jacobian.)
   */
  VectorD2 iterate(const VectorD2 &x0, size_t step_limit = 25)
  {
    steps = 0;
    ssc_steps = 0;
    VectorD2 x = x0;
    MatrixD2 J(dim, dim);

    for (size_t k = 0; k < step_limit; k++)
      {
        J = f.diff(x);    // Jacobian J
        J.gauss_jordan(); // Inverse J^{-1}
        step(J, x);       // Newton step

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

    throw std::invalid_argument("step limit exceeded");
  }

  /*! \fn iterate_broyden
   *  \brief Use rank-1 updates during iteration.
   *
   * The Jacobian is computed periodicially, as specified trough the \c skips
   * parameter. This method converges slower, with the default step limit is
   * chosen accordingly.
   *
   * For \b small \b problems such as the Thomas-Fermi problem, the slow convergence
   * rate was more expensive than computing the Jacobian in each step.
   *
   * For \b large \b problems, the Sherman-Morrison formula may be used to update
   * \f$J^{-1}\f$ directly, instead of performing an LU decomposition. Due to
   * lower relevance for small problems, this is not implemented here.
   */
  VectorD2 iterate_broyden(const VectorD2 &x0, size_t skips = 5,
                           size_t step_limit = 50)
  {
    VectorD2 x = x0;
    MatrixD2 J = f.diff(x);
    MatrixD2 J_inv(J);
    J_inv.gauss_jordan();

    steps = 1;
    ssc_steps = 0;

    // Set new value of x, preserve J_inv
    VectorD2 x_prev = x; //copy
    step(J_inv, x);

    for (size_t k = 1; k < step_limit; k++)
      {
        if (k % skips == 0)
          {
            // Compute exact Jacobian
            J = f.diff(x);
            J_inv = J;
            J_inv.gauss_jordan();

            x_prev = x; //copy
            step(J_inv, x);
          }
        else
          {
            VectorD2 p = x - x_prev;
            VectorD2 q = f(x) - f(x_prev);
            MatrixD2 S(p.size(), p.size());

            // Rank-1 updates: update J_n from J_{n-1}
            // MatrixD2 V(p.size(), p.size());
            // V.outer_product(q - J * p, p);
            // J = J + 1. / p.norm_sqr() * V;

            // Sherman-Morrison formula: update J_inv_n from J_inv_{n-1}
            S.outer_product((p - J_inv*q) / J_inv.matrix_scalar_product(p, q), p);
            J_inv = J_inv + S * J_inv;

            // Perform quasi-Newton step
            x_prev = x; //copy
            step(J_inv, x);
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

    throw std::invalid_argument("step limit exceeded");
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
