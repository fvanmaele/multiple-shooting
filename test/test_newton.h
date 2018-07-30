#ifndef TEST_NEWTON_H
#define TEST_NEWTON_H
#include <cmath>
#include <vector>

#include "../algo/newton.h"
#include "../base/forward_ad.h"

struct FAD_Nonlinear : public FAD_Functor
{
  virtual std::vector<FAD_Number>
  operator()(const std::vector<FAD_Number> &x) override
  {
    std::vector<FAD_Number> y = {
      std::pow(x[0],2) + std::pow(x[1],2) - 1,
      std::pow(x[0],2) - x[1]
    };
    return y;
  }
};

// Ortega, Scientific Computing, Tab. 4.3.2
void Test_Newton()
{
  FAD_Nonlinear f;
  Function_AD F(f, 2, 2);
  Newton N(F, 2);

  // Starting value
  dealii::Vector<FP_Type> s(2);
  s[0] = 0.5;
  s[1] = 0.5;

  size_t steps = 0;
  for (size_t k = 0; k < 20; k++)
    {
      s = N.step(F.jacobian(s), s, true);
      steps++;

      if (N.stopping_criterion(1e-8))
        {
          std::cout << "Solution of F: (" << steps << " steps) " << s;
          break;
        }
    }
}

#endif // TEST_NEWTON_H
