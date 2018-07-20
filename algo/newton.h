#ifndef NEWTON_H
#define NEWTON_H
#include "jacobian.h"
#include "../base/number_type.h"

class Newton_SSC
{
public:
  Newton_SSC(std::vector<FAD_Number> _f, std::vector<FAD_Number> _x0,
             size_t _nsteps = 1000) :
    f(_f), x(_x0), nsteps(_nsteps)
  {
    // By definition, f: R^d -> R^d
    assert(f.size() == x.size());
  }

  void iterate()
  {
    dealii::Vector<NumberType> d(f.size());
    for (size_i i = 0; i <= nsteps; i++)
      {
        FAD_Jacobian Nabla;
        dealii::FullMatrix<NumberType> J = Nabla(f, x.size());
      }
  }

private:
  std::vector<FAD_Number> f, x;
  size_t nsteps;
};

#endif // NEWTON_H
