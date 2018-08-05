#include <iostream>

#include "base/forward_ad.h"
#include "base/types.h"
#include "ivp/runge_kutta.h"
#include "test/run_tests.h"

VectorD2 ThomasFermi(FP_Type t, const VectorD2 &u)
{
  VectorD2 y(2);
  y[0] = t * u[1];
  y[1] = 4 * std::pow(u[0], 1.5);

  return y;
}

int main(int argc, char* argv[])
{
  all_tests();
  return 0;
}
