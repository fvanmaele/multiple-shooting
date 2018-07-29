#include <iostream>

#include "algo/newton.h"
#include "base/functor.h"
#include "base/forward_ad.h"
#include "base/types.h"
#include "ivp/runge_kutta.h"

#include "test/sheet1.h"
#include "test/sheet2.h"
#include "test/sheet3.h"
#include "test/test_forward_ad.h"
#include "test/test_newton.h"
#include "test/test_runge_kutta.h"

class ThomasFermi : public RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    dealii::Vector<FP_Type> result(2);
    result[0] = t * u[1];
    result[1] = 4 * std::pow(u[0], 1.5);

    return result;
  }
};

int main(int argc, char* argv[])
{
  Test_Sheet1();
  Test_Sheet2();
  Test_Sheet3();
  Test_FAD();
  std::cout << std::endl;
  Test_Newton();
  return 0;
}
