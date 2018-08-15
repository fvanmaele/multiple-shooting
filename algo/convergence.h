#ifndef CONVERGENCE_H
#define CONVERGENCE_H
#include <algorithm>
#include <cmath>
#include <vector>

#include "../algo/newton.h"
#include "../ivp/eos_method.h"
#include "../base/types.h"

class ConvergenceTester
{
public:
  ConvergenceTester(OneStepMethod &_M, size_t _max = 3, FP_Type _e = 1e-1) :
    M(_M), ROA(_max), TOL(_e)
  {}
  
  FP_Type eoc(const FP_Type t_lim, const FP_Type h)
  {
    std::vector<FP_Type> v;

    // A method is in its "asymptotic region of accuracy" when h is small enough
    // to give a good estimate of p. This required size of h can be different for
    // different problems. We make an estimate of p for several different h, and
    // check that we get approximately the same value.
    for (size_t i = 0; i < ROA; i++)
      {
	std::vector<dealii::Vector<FP_Type> > y;
	FP_Type h0 = h / std::pow(2., i);
	FP_Type h1 = h0 / 2.;
	FP_Type h2 = h0 / 4.;

	M.iterate_up_to(t_lim, h0);
	y.push_back(M.approx());

	M.iterate_up_to(t_lim, h1);
	y.push_back(M.approx());

	M.iterate_up_to(t_lim, h2);
	y.push_back(M.approx());

	// Difference between approximates of different step width
	const FP_Type norm1 = (y[0] - y[1]).l2_norm();
	const FP_Type norm2 = (y[1] - y[2]).l2_norm();

	FP_Type alpha = oc_alpha(norm1, norm2);
	v.push_back(alpha);
      }

    check_region(v, "EOC");
    return v.back();
  }

  FP_Type ooc(FP_Type t_lim, FP_Type h, const dealii::Vector<FP_Type> &u)                  
  {
    std::vector<FP_Type> v;

    for (size_t i = 0; i < ROA; i++)
      {
	std::vector<dealii::Vector<FP_Type> > y;
	FP_Type h0 = h / std::pow(2., i);
	FP_Type h1 = h0 / 2.;

	M.iterate_up_to(t_lim, h0);
	y.push_back(M.approx());

	M.iterate_up_to(t_lim, h1);
	y.push_back(M.approx());

	// Difference between approximates and exact value
	const FP_Type norm1 = (y[0] - u).l2_norm();
	const FP_Type norm2 = (y[1] - u).l2_norm();

	FP_Type alpha = oc_alpha(norm1, norm2);
	v.push_back(alpha);
      }

    check_region(v, "OOC");
    return v.back();
  }

private:
  OneStepMethod &M;
  size_t  ROA;
  FP_Type TOL;

  FP_Type oc_alpha(const FP_Type norm1, const FP_Type norm2)
  {
    return 1./std::log(2.) * std::log(norm1 / norm2);
  }

  FP_Type diameter(const std::vector<FP_Type> &v)
  {
    auto v_min = std::min_element(v.begin(), v.end());
    auto v_max = std::max_element(v.begin(), v.end());

    return std::fabs(*v_min - *v_max);
  }

  void check_region(const std::vector<FP_Type> &v, const char* str)
  {
    // Check if computed values are approximately equal by looking at the longest
    // distance between elements.
    FP_Type diam = diameter(v);

    if (diam > TOL)
      std::cerr << "warning: method not in asymptotic region of accuracy"
		<< std::endl;
    else
      std::cout << str << " steps: " << ROA
		<< std::endl
		<< str << " diameter: " << diam
		<< std::endl;
  }
};

#endif // CONVERGENCE_H
