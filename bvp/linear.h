#ifndef LINEAR_H
#define LINEAR_H

#include "../base/types.h"
#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

#include "shooting.h"

// Single shooting method for linear BVPs
//    y' = f(x,y),  A*y(a) + B*y(b) = c
//
// See Stoer, Num. Math. 2, pp.195
template <typename DiffMethod>
class LinearBVP : public DivFunctor
{
public:
  LinearBVP(TimeFunctor &_f, FP_Type _t0, FP_Type _t1,
            std::vector<FP_Type> _A,
            std::vector<FP_Type> _B,
            std::vector<FP_Type> _c) :
    f(_f), t0(_t0), t1(_t1), dim(_c.size())
  {
    assert(_A.size() == _B.size());
    assert(_A.size() == dim * dim);

    A = dealii::FullMatrix<FP_Type>(dim, dim, _A.data());
    B = dealii::FullMatrix<FP_Type>(dim, dim, _B.data());
    c = dealii::Vector<FP_Type>(_c.begin(), _c.end());
  }

  virtual dealii::Vector<FP_Type>
  operator()(const dealii::Vector<FP_Type> &s) override
  {
    DiffMethod S(f);
    dealii::Vector<FP_Type> y = S.solve_y(t0, t1, s);

    return A*s + B*y - c;
  }

  virtual dealii::FullMatrix<FP_Type>
  diff(const dealii::Vector<FP_Type> &s) override
  {
    DiffMethod S(f);
    dealii::FullMatrix<FP_Type> Z = S.solve_Z(t0, t1, s);

    return A + B*Z;
  }

private:
  TimeFunctor &f;
  FP_Type t0, t1;
  size_t dim;

  dealii::FullMatrix<FP_Type> A, B;
  dealii::Vector<FP_Type> c;
};

#endif // LINEAR_H
