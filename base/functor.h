#ifndef FUNCTOR_H
#define FUNCTOR_H
#include "../base/number_type.h"
#include "../lac/vector_operators.h"

class Functor
{
  // This class does not distinguish between scalars and vectors; in particular,
  // a scalar is represented as a vector of size 1.
public:
  virtual dealii::Vector<NumberType> operator()(NumberType t, const dealii::Vector<NumberType> &u) = 0;
  virtual ~Functor() = default;
};

// This functor represents the RHS of the autonomous IVP
//    u'(t) = u(t);
class RHS_P11 : public Functor
{
public:
  virtual dealii::Vector<NumberType> operator()(NumberType, const dealii::Vector<NumberType> &u) override
  {
    return u;
  }
};

// This functor represents the RHS of the IVP
//    u'(t) = t*u(t)
class RHS_P13 : public Functor
{
public:
  virtual dealii::Vector<NumberType> operator()(NumberType t, const dealii::Vector<NumberType> &u) override
  {
    return t*u;
  }
};

// This functor models the Lotka-Volterra equations.
class RHS_P21 : public Functor
{
public:
  // Constructor for the 4 parameters a, b, c, d.
  RHS_P21(NumberType _a, NumberType _b, NumberType _c, NumberType _d) :
    a(_a), b(_b), c(_c), d(_d)
  {}

  // u represents the vector u = (u1(t), u2(t)).
  virtual dealii::Vector<NumberType> operator()(NumberType t, const dealii::Vector<NumberType> &u) override
  {
    dealii::Vector<NumberType> result(2);
    result[0] = u[0] * (a - b * u[1]);
    result[1] = u[1] * (c * u[0] - d);
    return result;
  }

private:
  NumberType a, b, c, d;
};

#endif // FUNCTOR_H
