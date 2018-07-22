#ifndef RHS_H
#define RHS_H
#include "types.h"
#include "../lac/vector_operators.h"

// Note that there is no distinction between scalar and vector
// functions; both are represented as a dealii vector (of size 1
// in the scalar case).
class RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) = 0;
};

// This functor represents the RHS of the autonomous IVP
//    u'(t) = u(t) = f;
class RHS_P11 : public RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    return u;
  }
};

// This functor represents the RHS of the IVP
//    u'(t) = t*u(t) = f;
class RHS_P13 : public RHS
{
public:
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    return t * u;
  }
};

// This functor models the Lotka-Volterra equations.
class RHS_P21 : public RHS
{
public:
  // Constructor for the 4 parameters a, b, c, d.
  RHS_P21(FP_Type _a, FP_Type _b, FP_Type _c, FP_Type _d) :
    a(_a), b(_b), c(_c), d(_d)
  {}

  // u represents the vector u = (u1(t), u2(t)).
  virtual dealii::Vector<FP_Type>
  operator()(FP_Type t, const dealii::Vector<FP_Type> &u) override
  {
    dealii::Vector<FP_Type> result(2);
    result[0] = u[0] * (a - b * u[1]);
    result[1] = u[1] * (c * u[0] - d);

    return result;
  }

private:
  FP_Type a, b, c, d;
};

#endif // RHS_H
