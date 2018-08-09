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

class ShootingFunction
{
public:
  // Implementations for solve_y, solve_Z
  friend class SF_External;   // External differentation
  friend class SF_Automatic;  // Automatic differentation
  friend class SF_Manual;     // Use pre-computed fundamental matrix

  ShootingFunction(TimeFunctor &_f) : f(_f)
  {}

  // Compute y(t1; t0, s) by integrating IVP y(t0; s), where s represents
  // the specified initial value.
  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) = 0;

  // Compute D_s(y(t1; t0, s)) by external or exact differentation.
  // Return a pair of solutions, for solving Z typically involves solving y.
  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) = 0;

  virtual ~ShootingFunction() = default;

private:
  TimeFunctor &f;
};

class SF_External : public ShootingFunction
{
public:
  using ShootingFunction::ShootingFunction;

  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    AdaptiveMethod.iterate_with_ssc(t1, 1e-3, TOL, false);
    return AdaptiveMethod.approx();
  }

  // For the choice of TOL in the adaptive method and constant Epsilon,
  // see Stoer, Num. Math. 2, pp.192.
  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    dealii::FullMatrix<FP_Type> Z(s.size(), s.size());
    VectorD2 y = solve_y(t0, t1, s);

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
};

// Assumes f is DivFunctor with AD support.
class SF_Automatic : public ShootingFunction
{
public:
  using ShootingFunction::ShootingFunction;

  virtual VectorD2
  solve_y(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    // Compute solution of IVP
    AdaptiveMethod.iterate_with_ssc(t1, 1e-3, TOL, false);
    return AdaptiveMethod.approx();
  }

  virtual std::pair<VectorD2, MatrixD2>
  solve_Z(FP_Type t0, FP_Type t1, const VectorD2 &s) override
  {
    TimeDivFunctor* f_ad = dynamic_cast<TimeDivFunctor*>(&f);
    if (f_ad == nullptr)
      throw std::invalid_argument("functor is not differentiable");

    ERK<DOPRI87> AdaptiveMethod(f, t0, s);
    FP_Type TOL = std::sqrt(std::numeric_limits<FP_Type>::epsilon());

    // Compute IVP and variational equation simultaneously.
    // Step-size is controlled by the IVP only.
    AdaptiveMethod.iterate_with_ssc(t1, 1e-3, TOL, true);
    VectorD2 y = AdaptiveMethod.approx();
    MatrixD2 Z = AdaptiveMethod.fund_matrix();

    return std::make_pair(y, Z);
  }

private:
  // ShootingFunction::f
  // ShootingFunction::t0
  // ShootingFunction::t1
};

#endif // SHOOTING_H
