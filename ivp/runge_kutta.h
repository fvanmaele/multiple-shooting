#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>
#include <utility>

#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"
#include "eos_method.h"
#include "tableau.h"

/*!
 * \brief Solve an IVP of shape \f$u'(t) = f(t, u(t)), u(t_0) = u_0\f$
 * using an explicit Runge-Kutta method.
 */
template <typename ButcherTableau>
class ERK : public OneStepMethod
{
public:
  /*!
   * \brief ERK
   *
   * This constructor follows \c OneStepMethod. The used Butcher tableau is
   * specified via a template parameter. Whether a method is embedded (for
   * step-size control) is determined by checking if lower-order weights
   * are available in the Butcher tableau.
   *
   * Only fixed order methods or methods of order \f$p+1, p\f$ are supported.
   * See \c iterate_with_ssc for details.
   */
  ERK(TimeFunctor &f, FP_Type t0, VectorD2 u0, bool var_eq = false, Curve *u = nullptr) :
    OneStepMethod(f, t0, u0, var_eq, u)
  {
    ButcherTableau BTab;

    A  = MatrixD2(BTab.n, BTab.n, BTab.A.data());
    b1 = VectorD2(BTab.b_high.begin(), BTab.b_high.end());
    c  = VectorD2(BTab.c.begin(), BTab.c.end());
    p = BTab.p;

    if (BTab.b_low.size())
      {
        embedded_method = true;
        b2 = VectorD2(BTab.b_low.begin(), BTab.b_low.end());
        assert(b1.size() == b2.size());
      }
    else
      {
        embedded_method = false;
        b2 = VectorD2();
      }
  }

  VectorD2 k_increment(FP_Type t, const VectorD2 &y,
                       FP_Type h, const VectorD2 &b)
  {
    size_t s = b.size();
    size_t n = y.size();
    std::vector<VectorD2> k(s);
    k.at(0) = f(t, y);

    for (size_t i = 1; i < s; i++)
      {
        VectorD2 g(y.size());

        for (size_t j = 0; j < i; j++)
          g += A(i, j) * k.at(j);

        k.at(i) = f(t + h*c[i], y + h*g);
      }

    VectorD2 u_inc(n);

    for (size_t i = 0; i < s; i++)
      u_inc += b(i) * k.at(i);

    return u_inc;
  }


  std::pair<VectorD2, MatrixD2>
  k_variational(FP_Type t, const VectorD2 &y,
                FP_Type h, const VectorD2 &b, const MatrixD2 &Y)
  {
    size_t s = b.size();
    size_t n = y.size();
    assert(Y.m() == n);
    assert(Y.n() == n);

    // Increments k_1, ..., k_s for solution y(t)
    std::vector<VectorD2> k(s);
    k.at(0) = (*f_div)(t, y);

    // Increments K_1, ..., K_s for solution Y(t)
    std::vector<MatrixD2> K(s);
    K.at(0) = (*f_div).diff(t, y) * Y;

    for (size_t i = 1; i < s; i++)
      {
        VectorD2 g(n);
        MatrixD2 G(n);

        for (size_t j = 0; j < i; j++)
          {
            g += A(i, j) * k.at(j);
            G.add(1, A(i, j) * K.at(j));
          }

        k.at(i) = (*f_div)(t + h*c[i], y + h*g);
        K.at(i) = (*f_div).diff(t + h*c[i], y + h*g) * (Y + h*G);
      }

    VectorD2 inc_u(n);
    MatrixD2 inc_Y(n, n);

    for (size_t i = 0; i < s; i++)
      {
        inc_u += b(i) * k.at(i);
        inc_Y.add(1, b(i) * K.at(i));
      }
    return std::make_pair(inc_u, inc_Y);
  }

  virtual VectorD2
  increment_function(FP_Type t, const VectorD2 &y, FP_Type h) override
  { // Explicit method of higher order
    return k_increment(t, y, h, b1);
  }

  virtual std::pair<VectorD2, MatrixD2>
  increment_variational(FP_Type t, const VectorD2 &y, FP_Type h,
                        const MatrixD2 &Y) override
  {
    return k_variational(t, y, h, b1, Y);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  /*!
   * \brief iterate_with_ssc
   * \param t_lim
   * \param h0
   * \param TOL
   * \param C
   *
   * Iteration function with support for step-size control. In each step, two
   * solutions are computed, for order \f$p+1\f$ and \f$p\f$ respectively.
   * Whether a step is accepted is then decided by computing an "optimal"
   * step width. If the optimal width is smaller than the current value,
   * the step is repeated using the new width.
   *
   * More complex algorithms which choose method order dynamically (for example
   * the KARP method for orders 1 to 5) are not implemented, but may be more
   * suitable for discontinuous or rapidly changing ODEs. See the paper by Cash-Karp
   * for details.
   *
   * http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf
   *
   * When solving the variational equation, the step width is \e only determined
   * by the initial-value problem. This approach may be improved by first solving
   * the IVP with a \e continuous Runge-Kutta method. See the paper by L. Rández,
   * 1990.
   *
   * https://www.sciencedirect.com/science/article/pii/037704279290226N
   */
  void iterate_with_ssc(FP_Type t_lim, FP_Type h0, FP_Type TOL, FP_Type C = 2)
  {
    assert(embedded_method);
    size_t ad_count = 0;

    // Ensure initial step fits inside interval
    while (t0 + h0 > t_lim)
      {
        std::cerr << "warning: reducing initial step size (" << ad_count << ") ("
                  << h0 << "-> " << h0*1e-1 << ")" << std::endl;
        h0 *= 1e-1;
        ad_count++;

        if (ad_count > 3)
          throw std::invalid_argument("could not determine initial step size");
      }

    // Input variables
    FP_Type t = t0;
    FP_Type h_var = h0;
    size_t n = u0.size();

    // Output variables
    OneStepMethod::reset();
    misfires = 0;

    // Dynamic allocation, declare outside loop
    VectorD2 y = u0;
    VectorD2 inc_y1(n);
    VectorD2 inc_y2(n);
    MatrixD2 inc_Yn(n, n);

    // Algorithm 2.4.2
    while (t_lim - t > 0)
      { // (1) Candidates for u_k, v_k with time step h_k
        if (var_eq)
          {
            // The step-size is controlled only by the IVP, not the variational equation,
            // to avoid recomputing the Jacobian on a rejected time step.
            auto U = k_variational(t, y, h_var, b1, Yn);

            inc_y1 = y  + h_var * U.first;
            inc_y2 = y  + h_var * k_increment(t, y, h_var, b2);
            inc_Yn = Yn + h_var * U.second;
          }
        else
          {
            inc_y1 = y + h_var * k_increment(t, y, h_var, b1);
            inc_y2 = y + h_var * k_increment(t, y, h_var, b2);
          }
        steps++;

        if (sol_is_nan(inc_y1))
          throw std::overflow_error("global error too large (NaN)");

        // (2) Compute optimal step size.
        FP_Type local_error = (inc_y1 - inc_y2).l2_norm();
        FP_Type h_opt = h_var * std::pow(TOL / local_error, 1. / (p+1));
        //FP_Type h_opt = 0.9 * h_var * std::pow(TOL * std::abs(h_var) / local_error, 1./p);

        // (3) Time step is rejected; repeat step with optimal value.
        if (h_opt < h_var)
          {
            h_var = h_opt;
            misfires++;
            continue;
          }
        else
          { // (4) Time step was accepted.
            y  = inc_y1;
            Yn = inc_Yn;
            t += h_var;

            // Save intermediary solutions to IVP. (In the variational equation,
            // only the last step Y_n is stored.)
            OneStepMethod::save_step(t, y);
            step_sizes.push_back(h_var);

            // Set time step for next iteration.
            h_var = h_opt;

            // (*) Guarantee to hit right interval end (Rem. 2.4.3)
            if (t + 1.1*h_var >= t_lim)
              h_var = t_lim - t;

            // Comparison to exact solution
            if (u != nullptr)
              if (y.l2_norm() >= C*(*u)(t).l2_norm())
                throw std::domain_error("global error too large");
          }
      }

    if (timepoints.back() != t_lim)
      {
        std::string err = "time step outside interval end ("
            + std::to_string(timepoints.back()) + "; ["
            + std::to_string(t0) + ", "
            + std::to_string(t_lim) + "])";

        throw std::out_of_range(err.c_str());
      }
  }

  void print_step_size(std::ostream &out)
  {
    for (size_t i = 0; i < step_sizes.size(); i++)
      {
        out << timepoints[i] << "\t" << step_sizes[i] << std::endl;
      }
  }

private:
  // OneStepMethod::f;
  // OneStepMethod::u;
  // OneStepMethod::t0;
  // OneStepMethod::u0;
  // OneStepMethod::steps;
  // OneStepMethod::Yn
  // OneStepMethod::f_div;
  // OneStepMethod::timepoints;
  // OneStepMethod::uapprox;

  // Butcher tableau
  MatrixD2 A;
  VectorD2 b1, b2, c;
  size_t p;

  // Adaptive step-size plot
  std::vector<FP_Type> step_sizes;

  // Markers
  bool embedded_method;
  size_t misfires;
};

#endif // RUNGE_KUTTA_H
