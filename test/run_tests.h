#ifndef RUN_TESTS_H
#define RUN_TESTS_H

#include "sheet1.h"
#include "sheet2.h"
#include "sheet3.h"
#include "test_ad.h"
#include "test_bvp.h"
#include "test_newton.h"
#include "test_runge_kutta.h"
#include "test_vareq.h"

namespace Test
{
  void exercise_sheet()
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
  }

  void aut_diff()
  {
    Test_FAD();
    std::cout << std::endl
              << "== Newton Method =="
              << std::endl;
  }

  void linear_bvp()
  {
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

  void var_equation()
  {
    std::cout << "== Variational equation, problem 1 =="
              << std::endl;
    Test_P1();
    std::cout << std::endl
              << "== Variational equation, problem 2 =="
              << std::endl;
    Test_P2();
    std::cout << std::endl
              << "== Variational equation, problem 3 =="
              << std::endl;
    Test_P3();
    std::cout << std::endl
              << "== Variational equation, problem 4 =="
              << std::endl;
    Test_P4();
  }
}

#endif // RUN_TESTS_H
