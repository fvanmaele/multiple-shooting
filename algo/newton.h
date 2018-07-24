#ifndef NEWTON_H
#define NEWTON_H
#include "../base/types.h"
#include "../base/forward_ad.h"

template <size_t m, size_t n>
class Newton
{
public:
  Newton(FAD_Functor &f, const dealii::Vector<FP_Type> &x0)
    : F(f), x(x0), y(n), J(m, n)
  {
    // Nonlinear root finding problem with f: R^d -> R^d
    assert(x.size() == m);
    assert(y.size() == m);

    // Not optimized for systems of large dimension
    assert(x.size() <= 1000);
  }

  // Take the solution of the linear system Jd = y, and use it
  // to define the next step in the Newton method.
  void step(const dealii::Vector<FP_Type> &d,
            bool step_size_control = true,
            size_t ssc_limit = 20)
  {
    if (step_size_control)
      {
        // The norm (rhs) is fixed for each iteration j
        const FP_Type y_norm = y.norm_sqr();

        // Candidate for next step x_{k+1}
        dealii::Vector<FP_Type> x_cand(m);

        // Declare index outside of loop for bounds check
        size_t j;

        for (j = 0; j < ssc_limit; j++)
          {
            // Initialize new functor
            x_cand = x - pow(0.5, j)*d;
            F.init(x_cand);

            if (F.value().norm_sqr() < y_norm)
                break;
          }

        // j ist not guaranteed to be bounded; in practice, the step-size
        // control only triggers in the first few steps (see Remark 4.2.4)
        x = x_cand;
      }
    else
      {
        x -= d; // x_{k+1}
      }
  }

  void iterate(size_t steps = 1000,
               bool step_size_control = true,
               size_t ssc_limit = 20,
               bool rank1_updates = false)
  {
    if (rank1_updates)
      {
        assert(false); // to be done

        // Compute y_0, J_0 from starting value x_0
        x_prev = x;
        F.init(x_prev);

        y_prev = F.value();
        J_prev = F.jacobian();

        // Invert matrix for later use of Sherman-Morrison formula
        J_prev.invert();

        // Compute x_1 ... x_k, y_1 ... y_k, J_1 ... J_k
        for (size_t i = 2; i <= steps; i++)
          {
            //requires vector d (solution of linear system)
            //this->step(step_size_control);
            dealii::Vector<FP_Type> p = x - x_prev;
            dealii::Vector<FP_Type> q = y - y_prev;

            // ...
          }
      }
    else
      {
        for (size_t i = 1; i <= steps; i++)
          {
            F.init(x); // x_0 ... x_k
            xapprox.emplace_back(x);

            y = F.value();    // f (x_0) ... f (x_k)
            J = F.jacobian(); // f'(x_0) ... f'(x_k)

            // Solve linear system. We assume small dimensions d <= 1000,
            // and thus use a direct method (LU). The system is solved
            // in place; save a copy to compute j (equation 4.10) and
            // rank-1 updates (equation 4.14).
            dealii::Vector<FP_Type> d(y);
            J.compute_lu_factorization();
            J.solve(d);

            // Compute next iteration value x_{k+1}
            step(d, step_size_control, ssc_limit);
          }
      }
  }

  const dealii::Vector<FP_Type>
  approx() const
  {
    return x;
  }

  void reinit(const dealii::Vector<FP_Type> &x_new)
  {
    x = x_new;
    xapprox.clear();
  }

  void print(std::ostream &out)
  {
    for (size_t i = 0; i < xapprox.size(); i++)
      {
        out << i << "\t" << xapprox[i];
      }
  }

private:
  FAD_Wrapper<m, n> F;
  dealii::Vector<FP_Type> x;
  dealii::Vector<FP_Type> y;
  dealii::LAPACKFullMatrix<FP_Type> J;

  // Store previous step for rank-1 updates
  dealii::Vector<FP_Type> x_prev;
  dealii::Vector<FP_Type> y_prev;
  dealii::LAPACKFullMatrix<FP_Type> J_prev;

  // Save intermediary results
  std::vector<dealii::Vector<FP_Type> > xapprox;
};

#endif // NEWTON_H
