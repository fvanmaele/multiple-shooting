#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

class BoundaryCondition
{
public:
  virtual VectorD2
  operator()(const VectorD2 &u, const VectorD2 &v) = 0;

  virtual MatrixD2
  diff_u(const VectorD2 &u, const VectorD2 &v) = 0;

  virtual MatrixD2
  diff_v(const VectorD2 &u, const VectorD2 &v) = 0;

  virtual ~BoundaryCondition() = default;
};

class BC_Linear : public BoundaryCondition
{
public:
  BC_Linear(MatrixD2 _A, MatrixD2 _B, VectorD2 _c) :
    A(_A), B(_B), c(_c), dim(_c.size())
  {
    assert(A.n() == B.n());
    assert(A.m() == B.m());
    assert(A.n_elements() == dim*dim);
  }

  virtual VectorD2
  operator()(const VectorD2 &u, const VectorD2 &v) override
  {
    return A*u + B*v - c;
  }

  virtual MatrixD2
  diff_u(const VectorD2&, const VectorD2&) override
  {
    return A;
  }

  virtual MatrixD2
  diff_v(const VectorD2&, const VectorD2&) override
  {
    return B;
  }

private:
  MatrixD2 A, B;
  VectorD2 c;
  size_t dim;
};

#endif // BOUNDARY_H
