#ifndef SHOOTING_H
#define SHOOTING_H

#include "../base/functor.h"
#include "../base/types.h"
#include "../ivp/runge_kutta.h"

// Single shooting method for linear BVPs
//    y' = f(x,y),  A*y(a) + B*y(b) = c

// Implement shooting method using external differentation
class ShootingFunction_ED : public DivFunctor
{
public:
  ShootingFunction_ED(RHS &_f, FP_Type _t0, FP_Type _t1,
                      dealii::FullMatrix<FP_Type> _A,
                      dealii::FullMatrix<FP_Type> _B,
                      dealii::Vector<FP_Type> _c) :
    f(_f), t0(_t0), t1(_t1), A(_A), B(_B), c(_c)
  {}

  // XXX Implement caching?
  virtual dealii::Vector<FP_Type>
  value(const dealii::Vector<FP_Type> &s) override
  {
    // Difference to boundary value condition
    return A*s + B*eval_y(s) - c;
  }

  // For the choice of TOL in the adaptive method and constant Epsilon,
  // see Num. Math. 2, pp.192.
  virtual dealii::FullMatrix<FP_Type>
  jacobian(const dealii::Vector<FP_Type> &s) override
  {
    // DF(s) = A + B*Z(b;s)
    dealii::FullMatrix<FP_Type> Z(s.size(), s.size());

    for (size_t j = 0; j < s.size(); j++)
      {
        // Choice of delta
        FP_Type delta = std::sqrt(std::numeric_limits<FP_Type>::epsilon()) * s[j];

        // Value at (t1; s1 .. s_j+delta{s_j} .. s_n)
        dealii::Vector<FP_Type> s_e(s);
        s_e[j] += delta;

        // Value at (t1; s1 .. s_n)
        dealii::Vector<FP_Type> y = eval_y(s);

        // Approximation to partial derivative
        dealii::Vector<FP_Type> delta_y = (eval_y(s_e) - eval_y(s)) / delta;

        // Fill j-th column of "Jacobian" Z
        for (size_t i = 0; i < s.size(); i++)
          Z(i, j) = delta_y(i);
      }
    return A + B*Z;
  }

private:
  // Linear BVP
  RHS &f;
  FP_Type t0, t1;
  dealii::FullMatrix<FP_Type> A, B;
  dealii::Vector<FP_Type> c;

  // Value of y(t1; s) by integrating IVP y(t0; s)
  dealii::Vector<FP_Type>
  eval_y(const dealii::Vector<FP_Type> &s)
  {
    DOPRI Tab;
    ERK<7> AdaptiveMethod(f, t0, s, Tab.A, Tab.b1, Tab.b2, Tab.c);

    AdaptiveMethod.iterate_with_ssc(t1, 1e-1, 1e-8, 4);
    return AdaptiveMethod.approx();
  }
};

#endif // SHOOTING_H
