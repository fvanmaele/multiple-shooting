#ifndef SHOOTING_H
#define SHOOTING_H

#include <limits>
#include <utility>

#include "../algo/newton.h"
#include "../base/forward_ad.h"
#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"
#include "../ivp/runge_kutta.h"

/*!
 * \brief Integrate an IVP \f$u' = f(t, u)\f$ on a given time interval \f$[t_0, t_1]\f$
 * dependent on the initial value \f$s = u(t_0)\f$.
 *
 * \b Notation: \f$f(t, u(t; s))\f$ or \f$y(t; t_0, s)\f$.
 *
 * The partial derivatives \f$D_s f = \frac{\partial f}{\partial s}\f$ are computed
 * approximatively (\e external \e differentation) or by solving the variational
 * equation \f$Y' = \nabla_u f(t, u(t)) Y\f$.
 */
class ShootingFunction
{
public:
  // Implementations for solve_y, solve_Z
  template <typename M>
  friend class SF_External;   // External differentation
  template <typename M, typename N>
  friend class SF_Automatic;  // Automatic differentation

  /*!
   * \brief Constructor.
   *
   * When disabling step-size control, the intial step width \f$h_0\f$
   * should be set to a smaller value, for example \f$1e-3\f$.
   *
   * The appropriate value for \f$TOL\f$ depends on the chosen method
   * for differentiating \f$F\f$. For example, when computing \f$DF\f$
   * with external differentation, \f$F\f$ should be integrated as
   * accurately as possible.
   */
  ShootingFunction(TimeFunctor &_f, bool _ssc = true,
                   FP_Type _h0 = 1e-1, FP_Type _TOL = 1e-4)
    :
      f(_f), dim(_f.n_dim()), ssc(_ssc), h0(_h0), TOL(_TOL)
  {}

  /*!
   * \brief Return the dimension of \f$F(s)\f$.
   */
  size_t n_dim() const
  {
    return dim;
  }

  /*!
   * \brief Solve \f$y(t; t_0, s)\f$ in \f$t = t_1\f$.
   */
  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) = 0;

  /*!
   * \brief Solve \f$D_s y(t; t_0, s)\f$ in \f$t = t_1\f$.
   *
   * As solving \f$D_s\f$ typically involves solving \f$f\f$, a pair of
   * solutions is returned.
   */
  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) = 0;

  virtual ~ShootingFunction() = default;

private:
  TimeFunctor &f;
  size_t dim;
  bool ssc;
  FP_Type h0, TOL;
};

template <typename M>
class SF_External : public ShootingFunction
{
public:
  using ShootingFunction::ShootingFunction;

  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    ERK<M> Method(f, t0, s, false);

    if (ssc)
      Method.iterate_with_ssc(t1, h0, TOL);
    else
      Method.iterate_up_to(t1, h0);

    return Method.approx();
  }

  /*!
   * \brief Solve \f$D_s y(t; t_0, s)\f$ in \f$t = t_1\f$ by external differentation.
   *
   * For the choice of \f$TOL\f$ in the adaptive method and the constant
   * \f$eps\f$, see Stoer, Num. Math. 2, pp. 192.
   */
  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    dealii::FullMatrix<FP_Type> Z(s.size(), s.size());
    VectorD2 y = solve_y(t0, t1, s);

    for (size_t j = 0; j < s.size(); j++)
      { // Choice of delta
        FP_Type delta = std::sqrt(std::numeric_limits<FP_Type>::epsilon()) * s[j];

        // Use machine epsilon if j-th component is 0
        if (delta == 0)
          {
            std::cerr << "warning: s_j is zero" << std::endl;
            delta = std::numeric_limits<FP_Type>::epsilon();
          }

        // Value at (t1; s1 .. s_j+delta{s_j} .. s_n)
        VectorD2 s_e(s);
        s_e[j] += delta;

        // Approximation to partial derivative.
        VectorD2 delta_y = (solve_y(t0, t1, s_e) - y) / delta;

        // Fill j-th column of "Jacobian" Z.
        for (size_t i = 0; i < s.size(); i++)
          Z.set(i, j, delta_y(i));
      }

    return std::make_pair(y, Z);
  }

private:
  // ShootingFunction::f
  // ShootingFunction::ssc
  // ShootingFunction::h0
  // ShootingFunction::TOL
};

// Assumes f is DivFunctor with AD support.
template <typename M, typename M_Var = M>
class SF_Automatic : public ShootingFunction
{
public:
  using ShootingFunction::ShootingFunction;

  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    ERK<M> Method(f, t0, s, false);

    if (ssc)
      Method.iterate_with_ssc(t1, h0, TOL);
    else
      Method.iterate_up_to(t1, h0);

    return Method.approx();
  }

  /*!
   * \brief Solve \f$D_s y(t; t_0, s)\f$ in \f$t = t_1\f$ by automatic differentation.
   */
  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    TimeDivFunctor* f_ad = dynamic_cast<TimeDivFunctor*>(&f);

    if (f_ad == nullptr)
      throw std::invalid_argument("functor is not differentiable");

    // Compute IVP and variational equation simultaneously.
    ERK<M_Var> Method(f, t0, s, true);

    if (ssc)
      // Note: step-size is controlled by the IVP only.
      Method.iterate_with_ssc(t1, h0, TOL);
    else
      Method.iterate_up_to(t1, h0);

    VectorD2 y = Method.approx();
    MatrixD2 Z = Method.fund_matrix();
    return std::make_pair(y, Z);
  }

private:
  // ShootingFunction::f
  // ShootingFunction::ssc
  // ShootingFunction::h0
  // ShootingFunction::TOL
};

#endif // SHOOTING_H
