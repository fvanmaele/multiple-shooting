#include <iostream>

#include "algo/newton.h"
#include "base/forward_ad.h"
#include "base/types.h"
#include "ivp/runge_kutta.h"

#include "sheet/sheet1.h"
#include "sheet/sheet2.h"
#include "sheet/sheet3.h"
#include "test/test_forward_ad.h"
#include "test/test_newton.h"

int main(int argc, char* argv[])
{
//  Test_Sheet1();
//  Test_Sheet2();
//  Test_FAD();
//  std::cout << std::endl;
  Test_Newton();
//  Test_Sheet3();
  return 0;
}
