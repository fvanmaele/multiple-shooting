#ifndef RUN_TESTS_H
#define RUN_TESTS_H

#include "sheet1.h"
#include "sheet2.h"
#include "sheet3.h"
#include "test_ad.h"
#include "test_newton.h"
#include "test_runge_kutta.h"
#include "test_bvp.h"

void all_tests()
{
  std::cout << "== Exercise Sheet 1 =="
            << std::endl;
  Test_Sheet1();
  std::cout << std::endl
            << "== Exercise Sheet 2 =="
            << std::endl;
  Test_Sheet2();
  std::cout << std::endl
            << "== Exercise Sheet 3 =="
            << std::endl;
  Test_Sheet3();
  std::cout << std::endl
            << "== Automatic Differentation =="
            << std::endl;
  Test_FAD();
  std::cout << std::endl
            << "== Newton Method =="
            << std::endl;
  Test_Newton();
  std::cout << std::endl
            << "== BVP (Stoer, single shooting) =="
            << std::endl;
  Test_Stoer();
  std::cout << std::endl
            << "== BVP (Troesch, single shooting) =="
            << std::endl;
  Test_Troesch();
}

#endif // RUN_TESTS_H
