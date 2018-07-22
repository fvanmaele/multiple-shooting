#ifndef EOC_H
#define EOC_H
#include <vector>

#include "../base/functor.h"
#include "../base/types.h"

template <typename M>
FP_Type eoc(RHS &rhs, FP_Type t0, FP_Type t1, FP_Type h,
           dealii::Vector<FP_Type> u0)
{
  M Method1(rhs, t0, u0, h);
  M Method2(rhs, t0, u0, h/2);
  M Method3(rhs, t0, u0, h/4);

  // execute method
  Method1.iterate_up_to(t1);
  Method2.iterate_up_to(t1);
  Method3.iterate_up_to(t1);

  // retrieve difference between approximates of different step width
  dealii::Vector<FP_Type> diff1 = Method1.uapprox.back() - Method2.uapprox.back();
  dealii::Vector<FP_Type> diff2 = Method2.uapprox.back() - Method3.uapprox.back();

  // estimate EOC with given formula
  FP_Type Alpha = 1 / std::log(2.) * std::log(diff1.l2_norm() / diff2.l2_norm());
  return Alpha;
}

// This function take an IVP with known exact solution u, computes an approximate
// solution y (through a method M) at the point t1, and prints (the norm of) the
// difference between y(t1) and u(t1).
template <typename M>
void evaluate_with_eoc(RHS &rhs, FP_Type t0, FP_Type t1, FP_Type h,
                       const dealii::Vector<FP_Type> &u0,
                       const dealii::Vector<FP_Type> &u1)
{
  // Set up IVP method
  M Method(rhs, t0, u0, h);
  Method.iterate_up_to(t1);

  // approximate solution of IVP at t1
  dealii::Vector<FP_Type> y1 = Method.uapprox.back();
  size_t nsteps = Method.steps();

  // compute order of convergence
  FP_Type EOC = eoc<M>(rhs, t0, t1, h, u0);

  // accuracy of approximation with step length h
  FP_Type diff = (y1-u1).l2_norm();

  std::cout << "Evaluations: " << nsteps << std::endl
            << "EOC (h = " << h << "): " << EOC << std::endl
            << "Norm: " << diff << std::endl;
}

#endif // EOC_H
