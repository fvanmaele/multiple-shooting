#ifndef SHOOTING_H
#define SHOOTING_H

#include <limits>

#include "../algo/newton.h"
#include "../base/forward_ad.h"
#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"
#include "../ivp/runge_kutta.h"

// Single shooting method for linear BVPs
//    y' = f(x,y),  A*y(a) + B*y(b) = c
//
// See Stoer, Num. Math. 2, pp.195
class ShootingFunction : public DivFunctor
{
public:
  // Implementations for solve_y, solve_Z
  friend class SF_External;   // External differentation
  friend class SF_Automatic;  // Automatic differentation
  friend class SF_Manual;     // Use pre-computed fundamental matrix

  ShootingFunction(TimeFunctor &_f, FP_Type _t0, FP_Type _t1,
                   dealii::FullMatrix<FP_Type> _A,
                   dealii::FullMatrix<FP_Type> _B,
                   dealii::Vector<FP_Type> _c) :
    f(_f), t0(_t0), t1(_t1), A(_A), B(_B), c(_c)
  {}

  // A. Compute y(t1; s) by integrating IVP y(t0; s).
  //  * Called multiple times per Newton step through step-size control.
  virtual dealii::Vector<FP_Type>
  solve_y(const dealii::Vector<FP_Type> &s) = 0;

  // B. Compute d_y(t1; s)/d_s by external or exact differentation.
  //  * Requires solving the variational equation in each integration step,
  //    or approximating through difference quotients per component of s.
  //  * Called once per Newton step.
  virtual dealii::FullMatrix<FP_Type>
  solve_Z(const dealii::Vector<FP_Type> &s) = 0;

  virtual ~ShootingFunction() = default;

  virtual dealii::Vector<FP_Type>
  operator()(const dealii::Vector<FP_Type> &s) override
  {
    dealii::Vector<FP_Type> y = solve_y(s);
    return A*s + B*y - c;
  }

  virtual dealii::FullMatrix<FP_Type>
  diff(const dealii::Vector<FP_Type> &s) override
  {
    dealii::FullMatrix<FP_Type> Z = solve_Z(s);
    return A + B*Z;
  }

private:
  TimeFunctor &f;
  FP_Type t0, t1;
  dealii::FullMatrix<FP_Type> A, B;
  dealii::Vector<FP_Type> c;
};

class SF_External : public ShootingFunction
{
public:
  using ShootingFunction::ShootingFunction;

  virtual dealii::Vector<FP_Type>
  solve_y(const dealii::Vector<FP_Type> &s) override
  {
    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    AdaptiveMethod.iterate_with_ssc(t1, 1e-1, TOL, 7);
    return AdaptiveMethod.approx();
  }

  // For the choice of TOL in the adaptive method and constant Epsilon,
  // see Stoer, Num. Math. 2, pp.192.
  virtual dealii::FullMatrix<FP_Type>
  solve_Z(const dealii::Vector<FP_Type> &s) override
  {
    dealii::FullMatrix<FP_Type> Z(s.size(), s.size());
    dealii::Vector<FP_Type> y = solve_y(s);

    for (size_t j = 0; j < s.size(); j++)
      {
        // Choice of delta
        FP_Type delta = std::sqrt(std::numeric_limits<FP_Type>::epsilon()) * s[j];

        // Use machine epsilon if j-th component is 0
        if (delta == 0)
          {
            std::cerr << "warning: s_j is zero" << std::endl;
            delta = std::numeric_limits<FP_Type>::epsilon();
          }

        // Value at (t1; s1 .. s_j+delta{s_j} .. s_n)
        dealii::Vector<FP_Type> s_e(s);
        s_e[j] += delta;

        // Approximation to partial derivative.
        dealii::Vector<FP_Type> delta_y = (solve_y(s_e) - y) / delta;

        // Fill j-th column of "Jacobian" Z.
        for (size_t i = 0; i < s.size(); i++)
          Z.set(i, j, delta_y(i));
      }
    return Z;
  }

private:
  // ShootingFunction::f
  // ShootingFunction::t0
  // ShootingFunction::t1
  // ShootingFunction::A
  // ShootingFunction::B
  // ShootingFunction::c
};

// Assumes f is DivFunctor with AD support.
class SF_Automatic : public ShootingFunction
{
  using ShootingFunction::ShootingFunction;

  virtual dealii::Vector<FP_Type>
  solve_y(const dealii::Vector<FP_Type> &s) override
  {
    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    // Only compute solution of IVP
    AdaptiveMethod.iterate_with_ssc(t1, 1e-1, TOL, 7, false);
    return AdaptiveMethod.approx();
  }

  virtual dealii::FullMatrix<FP_Type>
  solve_Z(const dealii::Vector<FP_Type> &s) override
  {
    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    // Compute fundamental matrix simultaneously to solution of IVP
    AdaptiveMethod.iterate_with_ssc(t1, 1e-1, TOL, 7, true);
    return AdaptiveMethod.fund_matrix();
  }

private:
  // ShootingFunction::f
  // ShootingFunction::t0
  // ShootingFunction::t1
  // ShootingFunction::A
  // ShootingFunction::B
  // ShootingFunction::c
};

#endif // SHOOTING_H
