#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>

#include "../base/forward_ad.h"

// f: R^3 -> R^3,  n = m = 3
struct FAD_Test : public FAD_Functor
{
  virtual std::vector<FAD_Number>
  operator()(const std::vector<FAD_Number> &x) override
  {
    std::vector<FAD_Number> y = {
      std::pow(x[0], 3) * std::sin(x[1]),
      std::cos(x[0]) * std::sin(x[1]),
      std::exp(x[2])
    };

    return y;
  }
};

void Test_FAD()
{
  dealii::Vector<FP_Type> x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;

  // Define wrapper object
  FAD_Test f;
  FAD_Wrapper<3, 3> F(f);

  // Initialize templates and evaluate function
  F.init(x);
  dealii::Vector<FP_Type> y = F.value();
  std::cout << y << std::endl;

  // Evaluate Jacobian
  F.jacobian().print_formatted(std::cout, 3, true, 0, "0");
}

#endif // TEST_FORWARD_AD_H
