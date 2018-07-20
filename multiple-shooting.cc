#include <deal.II/lac/full_matrix.h>
#include <iostream>
#include "base/functor.h"
#include "base/number_type.h"
#include "ivp/runge_kutta.h"
#include "algo/jacobian.h"

typedef Sacado::Fad::DFad<double> fad_double;

int main(int argc, char* argv[])
{
  Test_FAD();
  return 0;
}
