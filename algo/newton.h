#ifndef NEWTON_H
#define NEWTON_H
#include "jacobian.h"
#include "../base/types.h"
#include "../base/functor.h"

template <size_t dim>
class Newton_SSC
{
public:
  Newton_SSC(dealii::Tensor<1, dim, FAD_Number> _f,
             dealii::Tensor<1, dim, FAD_Number> _x0,
             size_t _nsteps = 1000) :
    f(_f), x(_x0), nsteps(_nsteps)
  {
    // By definition, f: R^d -> R^d
    assert(f.dimension == x.dimension);
  }

  void iterate()
  {
    unsigned int n = f.dimension;
    for (size_i i = 0; i <= nsteps; i++)
      {
        FAD_Jacobian Nabla;
        dealii::LAPACKFullMatrix<FP_Type> J = Nabla(f, n);

        // Unroll values of f(x) and solve linear system
        dealii::Vector<FP_Type> d(n);
        for (size_t i = 0; i < n; i++)
          {
            b[i] = d[i].val();
          }

        J.compute_lu_factorization();
        J.solve(d);
        // Step size control

      }
  }

private:
  dealii::Tensor<1, dim, FAD_Number> f, x;
  size_t nsteps;
};

#endif // NEWTON_H
