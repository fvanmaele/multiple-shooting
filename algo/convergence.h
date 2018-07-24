#ifndef EOC_H
#define EOC_H
#include <algorithm>
#include <cmath>
#include <vector>

#include "../base/functor.h"
#include "../base/types.h"

template <typename M>
FP_Type eoc(RHS &rhs, FP_Type t0, FP_Type t1, FP_Type h,
            const dealii::Vector<FP_Type> &u0,
            size_t ROA_limit = 3, FP_Type TOL = 1e-1)
{
  std::vector<FP_Type> eoc_h;

  // A method is in its "asymptotic region of accuracy" when h is small enough
  // to give a good estimate of p. This required size of h can be different for
  // different problems. We make an estimate of p for several different h, and
  // check that we get approximately the same value.
  for (size_t i = 0; i < ROA_limit; i++)
    {
      M Method1(rhs, t0, u0, h / (1*std::pow(2, i)));
      M Method2(rhs, t0, u0, h / (2*std::pow(2, i)));
      M Method3(rhs, t0, u0, h / (4*std::pow(2, i)));

      // execute method
      Method1.iterate_up_to(t1);
      Method2.iterate_up_to(t1);
      Method3.iterate_up_to(t1);

      // retrieve difference between approximates of different step width
      dealii::Vector<FP_Type> Diff1 = Method1.approx() - Method2.approx();
      dealii::Vector<FP_Type> Diff2 = Method2.approx() - Method3.approx();

      FP_Type alpha = std::log(Diff1.l2_norm() / Diff2.l2_norm()) / std::log(2.);
      eoc_h.push_back(alpha);
    }

  // Check if computed values are approximately equal by looking at the longest
  // distance between elements.
  auto eoc_min = std::min_element(eoc_h.begin(), eoc_h.end());
  auto eoc_max = std::max_element(eoc_h.begin(), eoc_h.end());

  FP_Type diameter = std::fabs(*eoc_min - *eoc_max);
  std::cout << "EOC steps: " << ROA_limit  << std::endl
            << "EOC diameter: " << diameter << std::endl;

  if (diameter > TOL)
    throw std::domain_error("Could not determine EOC - is h small enough?");
  else
    return eoc_h.back();
}

template <typename M>
FP_Type ooc(RHS &rhs, FP_Type t0, FP_Type t1, FP_Type h,
            const dealii::Vector<FP_Type> &u0,
            const dealii::Vector<FP_Type> &u)
{
  M Method1(rhs, t0, u0, h);
  M Method2(rhs, t0, u0, h/2);

  Method1.iterate_up_to(t1);
  Method2.iterate_up_to(t1);

  // retrieve differences between approximates and exact solution
  dealii::Vector<FP_Type> Diff1 = Method1.approx() - u;
  dealii::Vector<FP_Type> Diff2 = Method2.approx() - u;

  return std::log(Diff1.l2_norm() / Diff2.l2_norm()) / std::log(2.);
}

// This function take an IVP with known exact solution u, computes an approximate
// solution y (through a method M) at the point t1, and prints (the norm of) the
// difference between y(t1) and u(t1).
template <typename M>
void evaluate_convergence(RHS &rhs, FP_Type t0, FP_Type t1, FP_Type h,
                          const dealii::Vector<FP_Type> &u0,
                          const dealii::Vector<FP_Type> &u)
{
  // Set up IVP method
  M Method(rhs, t0, u0, h);
  Method.iterate_up_to(t1);

  // approximate solution of IVP at t1
  dealii::Vector<FP_Type> y = Method.approx();
  size_t steps = Method.n_steps();

  // compute order of convergence
  FP_Type EOC = eoc<M>(rhs, t0, t1, h, u0);

  // accuracy of approximation with step length h
  FP_Type diff = (y - u).l2_norm();

  std::cout << "Evaluations: " << steps << std::endl
            << "EOC (h = " << h << "): " << EOC << std::endl
            << "Norm: " << diff << std::endl;
}

#endif // EOC_H
