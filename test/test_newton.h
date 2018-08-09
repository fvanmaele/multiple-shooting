#ifndef TEST_NEWTON_H
#define TEST_NEWTON_H
#include <cmath>
#include <vector>

#include "../algo/newton.h"
#include "../base/forward_ad.h"
#include "../lac/lac_types.h"

namespace Test
{
  VectorAD AutFAD(NumberAD t, const VectorAD &x)
  {
    VectorAD y = {
      std::pow(x[0], 2) + std::pow(x[1], 2) - 1,
      std::pow(x[0], 2) - x[1]
    };
    return y;
  }

  // Ortega, Scientific Computing, Tab. 4.3.2
  void Ortega()
  {
    // Starting value
    dealii::Vector<FP_Type> s(2);
    s[0] = 0.5;
    s[1] = 0.5;

    FAD_cWrapper F(AutFAD, 2);
    Newton N(F, 2);

    N.iterate(s);
  }
}

#endif // TEST_NEWTON_H
