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
class SingleShooting : public DivFunctor
{
public:
  SingleShooting(TimeFunctor &_f, FP_Type _a, FP_Type _b,
                 BoundaryCondition &_r) :
    f(_f), a(_a), b(_b), r(_r)
  {}

  // Called multiple times per Newton step through step-size control.
  virtual VectorD2
  operator()(const VectorD2 &s) override
  {
    DiffMethod M(f);
    VectorD2 y = M.solve_y(a, b, s);

    return r(s, y);
  }

  // Called once per Newton step.
  virtual MatrixD2
  diff(const VectorD2 &s) override
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
class MultipleShooting : public DivFunctor
{
public:
  MultipleShooting(TimeFunctor &_f, FP_Type _a, FP_Type _b,
                   BoundaryCondition &_r,
                   std::vector<FP_Type> _t,
                   std::vector<VectorD2> _s0) :
    f(_f), a(_a), b(_b), r(_r), t(_t), s0(_s0)
  {}

  virtual VectorD2
  operator()(const VectorD2 &s) override
  {
    size_t m = t.size();
    size_t n = s0.front().size();
    DiffMethod M(f);

    // Construct F(s)
    dealii::BlockVector<FP_Type> F(m, n); // F_1 ... F_m

    for (size_t i = 0; i < m; i++)
      if (i == m-1)
        F.block(i) = r(s.at(0), s.at(i));
      else
        F.block(i) = M.solve_y(t.at(i), t.at(i+1), s.at(i)) - s.at(i+1);
    return F;
  }

  virtual MatrixD2
  diff(const VectorD2 &s) override
  {
    size_t m = t.size();
    size_t n = s0.front().size();
    DiffMethod M(f);

    // XXX: use BlockMatrix or BlockSparseMatrix (Trilinos)
    MatrixD2 DF(m*n, m*n);
    MatrixD2 I = dealii::IdentityMatrix(n);

    for (size_t i = 0; i < m; i++)
      if (i == m-1)
        {
          DF.add(r.diff_u(s.at(0), s.at(i)), 1, i*n, 0);
          DF.add(r.diff_v(s.at(0), s.at(i)), 1, i*n, i*n);
        }
      else
        {
          std::pair<VectorD2, MatrixD2> G = M.solve_Z(t.at(i), t.at(i+1), s.at(i));
          MatrixD2 Z = G.second;

          DF.add(Z,  1, i*n, i*n);
          DF.add(I, -1, i*n, (i+1)*n);
        }
    return DF;
  }

private:
  TimeFunctor &f;           // rhs of ODE
  FP_Type a, b;             // interval boundaries
  BoundaryCondition &r;     // r(u, v)
  std::vector<FP_Type> t;   // interval subdivision
  std::vector<VectorD2> s0; // starting values
};

#endif // METHODS_H
