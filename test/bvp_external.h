#ifndef BVP_EXTERNAL_H
#define BVP_EXTERNAL_H
#include "../base/functor.h"
#include "../base/types.h"

// RHS for the initial value problem:
//    u' = (u1, u2)' = (u2, 1.5 * u1)
class RHS_BVP1 : public RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u)
  {
    assert(u.size() == 2);
    dealii::Vector<FP_Type> result(2);

    result[0] = u[1];
    result[1] = 1.5 * u[0];

    return result;
  }
};

// Stoer, Bulirsch, Num. Math 2, pp.192
void Test_SB192()
{
  RHS_BVP1 f;
  FP_Type t0 = 0.0;
  FP_Type t1 = 1.0;
}

#endif // BVP_EXTERNAL_H
