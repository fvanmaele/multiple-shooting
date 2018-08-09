#ifndef METHODS_H
#define METHODS_H

#include <deal.II/lac/block_vector.h>

#include "shooting.h"
#include "boundary.h"

// Single shooting method for BVPs
//    y' = f(x,y),  r(y(a), y(b)) = 0
//
// See Stoer, Num. Math. 2, pp.195
template <typename DiffMethod>
class SingleShooting
{
public:
  SingleShooting(TimeFunctor &_f, FP_Type _a, FP_Type _b,
                 BoundaryCondition &_r) :
    f(_f), a(_a), b(_b), r(_r)
  {}

  // Called multiple times per Newton step through step-size control.
  VectorD2 operator()(const VectorD2 &s)
  {
    DiffMethod M(f);
    VectorD2 y = M.solve_y(a, b, s);

    return r(s, y);
  }

  // Called once per Newton step.
  MatrixD2 diff(const VectorD2 &s)
  {
    DiffMethod M(f);
    std::pair<VectorD2, MatrixD2> G = M.solve_Z(a, b, s);

    VectorD2 y = G.first;
    MatrixD2 Z = G.second;

    return r.diff_u(s, y) + r.diff_v(s, y)*Z;
  }

private:
  TimeFunctor &f;
  FP_Type a, b;
  BoundaryCondition &r;
};

template <typename DiffMethod>
class MultipleShooting
{
public:
  MultipleShooting(TimeFunctor &_f, std::vector<FP_Type> _t,
                   BoundaryCondition &_r) :
    f(_f), t(_t), r(_r)
  {}

  dealii::BlockVector<FP_Type>
  operator()(const dealii::BlockVector<FP_Type> &s)
  {
    size_t m = t.size();
    // XXX: assumes all blocks have the same dimension
    size_t n = s.block(0).size();
    DiffMethod M(f);

    // Construct F(s)
    dealii::BlockVector<FP_Type> F(m, n); // F_1 ... F_m

    for (size_t i = 0; i < m; i++)
      if (i == m-1)
        F.block(i) = r(s.block(0), s.block(i));
      else
        F.block(i) = M.solve_y(t.at(i), t.at(i+1), s.block(i)) - s.block(i+1);

    return F;
  }

  dealii::FullMatrix<FP_Type>
  diff(const dealii::BlockVector<FP_Type> &s)
  {
    size_t m = t.size();
    // XXX: assumes all blocks have the same dimension
    size_t n = s.block(0).size();
    DiffMethod M(f);

    // XXX: use TrilinosWrappers::BlockSparseMatrix?
    MatrixD2 DF(m*n, m*n);
    MatrixD2 I = dealii::IdentityMatrix(n);

    for (size_t i = 0; i < m; i++)
      if (i == m-1)
        {
          DF.add(r.diff_u(s.block(0), s.block(i)), 1, i*n, 0);
          DF.add(r.diff_v(s.block(0), s.block(i)), 1, i*n, i*n);
        }
      else
        {
          std::pair<VectorD2, MatrixD2> G = M.solve_Z(t.at(i), t.at(i+1), s.block(i));
          MatrixD2 Z = G.second;

          DF.add(Z,  1, i*n, i*n);
          DF.add(I, -1, i*n, (i+1)*n);
        }

    return DF;
  }

private:
  TimeFunctor &f;           // rhs of ODE
  std::vector<FP_Type> t;   // interval subdivision
  BoundaryCondition &r;     // r(u, v)
};

#endif // METHODS_H
