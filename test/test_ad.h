#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>
#include <array>

#include "../base/forward_ad.h"

namespace Test
{
  // f: R^2 -> R^3,  m = 2, n = 3
  VectorAD FAD_23dim(NumberAD, const VectorAD &u)
  {
    assert(u.size() == 2);
    VectorAD y(3);

    y[0] = std::pow(u[0], 3) * std::sin(u[1]);
    y[1] = std::cos(u[0]) * std::sin(u[1]);
    y[2] = std::exp(u[1]);

    return y;
  }

  void FAD()
  {
    FP_Type t = 0;
    VectorD2 u(2);
    u[0] = 1;
    u[1] = 2;

    // Define wrapper object
    FAD_Setup F(FAD_23dim, 2, 3);

    // Initialize templates and evaluate function
    F.init(t, u);
    VectorD2 y = F.value();
    std::cout << y << std::endl;

    // Evaluate Jacobian
    F.diff().print_formatted(std::cout, 3, true, 0, "0");
  }
}

#endif // TEST_FORWARD_AD_H
