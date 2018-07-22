#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>

#include "../base/forward_ad.h"

typedef std::vector<FAD_Number> FAD_Vector;

// f: R^3 -> R^3,  n = m = 3
struct FAD_Test : public FAD_Functor
{
  virtual FAD_Vector
  operator()(const FAD_Vector &x) override
  {
    FAD_Vector y(3);

    y[0] = std::pow(x[0], 3) * std::sin(x[1]);
    y[1] = std::cos(x[0]) * std::sin(x[1]);
    y[2] = std::exp(x[2]);

    return y;
  }
};

void Test_FAD()
{
  dealii::Vector<FP_Type> x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;

  FAD_Test f;
  FAD_Wrapper<3, 3> F(f);
  F.init(x);

  std::cout << F.value() << std::endl;
  F.jacobian().print_formatted(std::cout, 3, true, 0, "0");
}

#endif // TEST_FORWARD_AD_H
