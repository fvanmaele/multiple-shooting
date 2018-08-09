#ifndef TEST_FORWARD_AD_H
#define TEST_FORWARD_AD_H
#include <cmath>
#include <array>

#include "../base/forward_ad.h"

namespace Test
{
  // f: R^3 -> R^3,  n = m = 3
  VectorAD FAD_3dim(NumberAD t, const VectorAD &u)
  {
    VectorAD y(3);
    y[0] = std::pow(u[0], 3) * std::sin(u[1]);
    y[1] = std::cos(u[0]) * std::sin(u[1]);
    y[2] = std::exp(u[2]);

    return y;
  }

  void FAD()
  {
    FP_Type t = 0;
    VectorD2 u(3);
    u[0] = 1;
    u[1] = 2;
    u[2] = 3;

    // Define wrapper object
    FAD_tWrapper F(FAD_3dim, 3);

    // Initialize templates and evaluate function
    VectorD2 y = F(t, u);
    std::cout << y << std::endl;

    // Evaluate Jacobian
    F.diff().print_formatted(std::cout, 3, true, 0, "0");
  }
}

#endif // TEST_FORWARD_AD_H
