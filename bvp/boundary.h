#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../lac/lac_types.h"
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"

/*! \class BoundaryCondition
 * \brief Abstract base class for boundary conditions \f$r(u,v)\f$
 * in boundary value problems.
 *
 * Supports evaluation and differentation over \f$u\f$ and \f$v\f$.
 */
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

/*! \class BC_Linear
 * \brief Boundary conditions for linear boundary value problems.
 *
 * A linear BVP \f$u' = f(t, u)\f$ has boundary condition \f$r(u,v) = Au + Bv - c\f$,
 * where \f$A\f$ and \f$B\f$ are quadratic, \f$n^2\f$-dimensional matrices.
 *
 * On the boundary of the interval \f$I = [a,b]\f$, there holds \f$r(u(a), u(b)) = 0\f$.
 * If \f$n = 2\f$, \f$A_{1}y(a) = c_1\f$ and \f$B_{2}y(b) = c_2\f$, we say the
 * BVP is \b separated.
 *
 * Most BVPs we consider are separated, for example the Thomas-Fermi or Troesch problems.
 */
class BC_Linear : public BoundaryCondition
{
public:
  /*! \fn BC_Linear
   * \brief Constructor.
   */
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
