#ifndef METHODS_H
#define METHODS_H

#include <algorithm> // std::iota

#include "shooting.h"
#include "boundary.h"

// Single shooting method for BVPs
//    y' = f(x,y),  r(y(a), y(b)) = 0
//
// See Stoer, Num. Math. 2, pp.195
class SingleShooting : public DivFunctor
{
public:
  SingleShooting(ShootingFunction &_M, FP_Type _a, FP_Type _b, BoundaryCondition &_r)
    : DivFunctor(_M.n_dim()), M(_M), a(_a), b(_b), r(_r)
  {}

  // Called multiple times per Newton step through step-size control.
  virtual VectorD2 operator()(const VectorD2 &s) override
  {
    VectorD2 y = M.solve_y(a, b, s);

    return r(s, y);
  }

  // Called once per Newton step.
  virtual MatrixD2 diff(const VectorD2 &s) override
  {
    auto G = M.solve_Z(a, b, s);
    VectorD2 y = G.first;
    MatrixD2 Z = G.second;

    return r.diff_u(s, y) + r.diff_v(s, y)*Z;
  }

private:
  ShootingFunction &M;
  FP_Type a, b;
  BoundaryCondition &r;
};

class MultipleShooting : public DivFunctor
{
public:
  MultipleShooting(ShootingFunction &_M, std::vector<FP_Type> _t, BoundaryCondition &_r)
    : DivFunctor(_M.n_dim()), M(_M), t(_t), r(_r), m(t.size())
  {}

  virtual VectorD2 operator()(const VectorD2 &s) override
  {
    size_t dim = n_dim();
    assert(s.size() == m * dim);
    VectorD2 F(s.size());

    for (size_t i = 0; i < m; i++)
      if (i == m - 1)
        { // Extract first block s_1 of s
          VectorD2 s_1(dim);
          for (size_t k = 0; k < dim; k++)
            s_1[k] = s[k];

          // Extract last block s_m of s
          VectorD2 s_m(dim);
          for (size_t k = 0; k < dim; k++)
            s_m[k] = s[k + (m-1)*dim];

          // Fill last row (boundary condition)
          VectorD2 r_i = r(s_1, s_m);
          for (size_t k = 0; k < dim; k++)
            F[k + (m-1)*dim] = r_i[k];
        }
      else
        { // Extract i-th block of s
          VectorD2 s_i(dim);
          for (size_t k = 0; k < dim; k++)
            s_i[k] = s[k + i*dim];

          // Extract (i+1)th block of s.
          VectorD2 s_j(dim);
          for (size_t k = 0; k < dim; k++)
            s_j[k] = s[k + (i+1)*dim];

          // Solve IVP and fill i-th component
          VectorD2 F_i = M.solve_y(t.at(i), t.at(i+1), s_i) - s_j;
          for (size_t k = 0; k < dim; k++)
            F[k + i*dim] = F_i[k];
        }

    return F;
  }

  virtual MatrixD2 diff(const VectorD2 &s) override
  {
    size_t dim = n_dim();
    assert(s.size() == m * dim);

    // XXX: use TrilinosWrappers::BlockSparseMatrix?
    MatrixD2 DF(s.size(), s.size());
    MatrixD2 I = dealii::IdentityMatrix(dim);

    for (size_t i = 0; i < m; i++)
      if (i == m - 1)
        { // Extract first block s_1 of s
          VectorD2 s_1(dim);
          for (size_t k = 0; k < dim; k++)
            s_1[k] = s[k];

          // Extract last block s_m of s
          VectorD2 s_m(dim);
          for (size_t k = 0; k < dim; k++)
            s_m[k] = s[k + (m-1)*dim];

          // Fill last row (boundary condition derivative)
          DF.add(r.diff_u(s_1, s_m), 1, i*dim, 0);
          DF.add(r.diff_v(s_1, s_m), 1, i*dim, i*dim);
        }
      else
        { // Extract i-th block of s
          VectorD2 s_i(dim);
          for (size_t k = 0; k < dim; k++)
            s_i[k] = s[k + i*dim];

          // Solve variational equation
          auto G = M.solve_Z(t.at(i), t.at(i+1), s_i);
          MatrixD2 Z = G.second;

          // Fill matrix
          DF.add(Z,  1, i*dim, i*dim);
          DF.add(I, -1, i*dim, (i+1)*dim);
        }

    return DF;
  }

private:
  ShootingFunction &M;
  std::vector<FP_Type> t;
  BoundaryCondition &r;
  size_t m;
};

#endif // METHODS_H
