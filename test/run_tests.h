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
    Test::Sheet1();
    std::cout << std::endl
              << "== Exercise Sheet 2 =="
              << std::endl;
    Test::Sheet2();
    std::cout << std::endl
              << "== Exercise Sheet 3 =="
              << std::endl;
    Test::Sheet3();
    std::cout << std::endl
              << "== Automatic Differentation =="
              << std::endl;
  }

  void aut_diff()
  {
    Test::FAD();
  }

  void linear_bvp()
  {
    std::cout << "== Newton Method =="
              << std::endl;
    Test::Ortega();
    std::cout << std::endl
              << "== BVP (Stoer, single shooting) =="
              << std::endl;
    Test::Stoer();
    std::cout << std::endl
              << "== BVP (Troesch, single shooting) =="
              << std::endl;
    Test::Troesch();
    std::cout << std::endl;
  }

  void var_equation()
  {
    std::cout << "== Variational equation, problem 1 =="
              << std::endl;
    Test::Var1();
    std::cout << std::endl
              << "== Variational equation, problem 2 =="
              << std::endl;
    Test::Var2();
    std::cout << std::endl
              << "== Variational equation, problem 3 =="
              << std::endl;
    Test::Var3();
    std::cout << std::endl
              << "== Variational equation, problem 4 =="
              << std::endl;
    Test::Var4();
  }
}

#endif // RUN_TESTS_H
