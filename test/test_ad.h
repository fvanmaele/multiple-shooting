#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>
#include <array>

#include "../base/forward_ad.h"

// f: R^3 -> R^3,  n = m = 3
struct FAD_Test : public FAD_TdVecField
{
  virtual std::vector<FAD_Number>
  operator()(FAD_Number t, const std::vector<FAD_Number> &u) override
  {
    std::vector<FAD_Number> y = {
      std::pow(u[0], 3) * std::sin(u[1]),
      std::cos(u[0]) * std::sin(u[1]),
      std::exp(u[2])
    };

    return y;
  }
};

void Test_FAD()
{
  FP_Type t = 1;
  dealii::Vector<FP_Type> u(3);
  u[0] = 1;
  u[1] = 2;
  u[2] = 3;

  // Define wrapper object
  FAD_Test f;
  FAD_Wrapper<3> F(f);

  // Initialize templates and evaluate function
  F.init(t, u);
  dealii::Vector<FP_Type> y = F.value();
  std::cout << y << std::endl;

  // Evaluate Jacobian
  F.nabla_u().print_formatted(std::cout, 3, true, 0, "0");
}

#endif // TEST_FORWARD_AD_H
