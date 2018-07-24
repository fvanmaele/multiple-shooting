#ifndef TEST_NEWTON_H
#define TEST_NEWTON_H

#include "../algo/newton.h"
#include "../base/forward_ad.h"
#include "../lac/vector_operators.h"

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

void Test_Newton()
{
  FAD_Nonlinear f;
  Newton<2,2> F(f, STV({0.5, 0.5}));

  F.iterate(20);
  F.print(std::cout);
}

#endif // TEST_NEWTON_H
