#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>
#include <array>

#include "../base/forward_ad.h"

namespace Test
{
  // f: R^3 -> R^3,  n = m = 3
  VectorAD FAD_3dim(const VectorAD &u)
  {
    VectorAD y(3);
    y[0] = std::pow(u[0], 3) * std::sin(u[1]);
    y[1] = std::cos(u[0]) * std::sin(u[1]);
    y[2] = std::exp(u[2]);

    return y;
  }

  void FAD()
  {
    VectorD2 u(3);
    u[0] = 1;
    u[1] = 2;
    u[2] = 3;

    // Define wrapper object
    FAD_cWrapper F(FAD_3dim, 3);

    // Initialize templates and evaluate function
    F.init(u);
    VectorD2 y = F.value();
    std::cout << y << std::endl;

    // Evaluate Jacobian
    F.diff().print_formatted(std::cout, 3, true, 0, "0");
  }
}

#endif // TEST_FORWARD_AD_H
