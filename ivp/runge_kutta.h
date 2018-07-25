#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>

#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"
#include "butcher_tableau.h"
#include "eos_method.h"

class ERK : public EOS_Method
{
public:
  // In addition to the usual arguments, the Range_Kuuta method
  // requires a matrix A and two vectors b, c. The matrix A is stored
  // as a dense matrix, for simplicity and due to its small dimension.
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::LAPACKFullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b,
      dealii::Vector<FP_Type> _c) :
    EOS_Method(f, t0, u0), A(_A), b1(_b), c(_c)
  {
    assert(c.size() == b1.size());
    embedded_method = false;
  }

  // Constructor for embedded methods
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::LAPACKFullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b1,
      dealii::Vector<FP_Type> _b2,
      dealii::Vector<FP_Type> _c) :
    EOS_Method(f, t0, u0), A(_A), b1(_b1), b2(_b2), c(_c)
  {
    assert(c.size() == b1.size());
    assert(c.size() == b2.size());
    embedded_method = true;
  }

  std::vector<dealii::Vector<FP_Type> >
  k_stage(const FP_Type &t, const dealii::Vector<FP_Type> &u,
          const FP_Type &h)
  {
    // Optional, see Remark 2.3.4
    std::vector<dealii::Vector<FP_Type> > g(c.size());
    // Keep values k per step for computation of y-hat
    std::vector<dealii::Vector<FP_Type> > k(c.size());

    // Compute: k_1 ... k_s
    for (size_t i = 0; i < c.size(); i++)
      {
        dealii::Vector<FP_Type> g_sum(u.size());

        // Compute: k_1 ... k_{i-1}
        // When i = 0, this loop is not executed, i.e. g_sum remains 0.
        for (size_t j = 0; j < i; j++)
          {
            FP_Type a_ij = A(i, j);

            // While the upper half of the matrix A is ignored,
            // entries in the lower half may still be zero.
            if (a_ij > 0)
              g_sum += a_ij * k.at(j);
          }

        g.at(i) = u + h*g_sum;
        k.at(i) = f(t + h*c(i), g.at(i));
      }
    return k;
  }

  dealii::Vector<FP_Type>
  k_increment(const FP_Type &t, const dealii::Vector<FP_Type> &u,
              const FP_Type &h, const dealii::Vector<FP_Type> &b)
  {
    // Compute stages k_1,..,k_s
    std::vector<dealii::Vector<FP_Type> > k = k_stage(t, u, h);
    // Initialize sum vector
    dealii::Vector<FP_Type> result(u.size());

    for (size_t i = 0; i < c.size(); i++)
      result += b(i) * k.at(i);
    return result;
  }

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    // Explicit method of higher order
    return k_increment(t, u, h, b1);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  void iterate_with_ssc(FP_Type t_lim, FP_Type h0 = 1e-1,
                        FP_Type TOL = 1e-4, size_t order = 5)
  {
    assert(embedded_method);
    FP_Type t = t0;     // start time
    FP_Type h_var = h0; // initial step size

    // Dynamic allocation, declare outside loop
    dealii::Vector<FP_Type> u = u0;
    dealii::Vector<FP_Type> inc_u(u.size());
    dealii::Vector<FP_Type> inc_v(u.size());

    // Init output variables
    EOS_Method::reset();
    steps = 0;
    misfires = 0;

    // Algorithm 2.4.2
    while (t_lim - t > 0)
      { // (1) Candidates for u_k, v_k with time step h_k
        inc_u = k_increment(t, u, h_var, b1);
        inc_v = k_increment(t, u, h_var, b2);
        steps++;

        // (2) Compute optimal step size.
        // Use norm properties to save operations.
        FP_Type local_error = std::abs(h_var) * (inc_u - inc_v).l2_norm();
        FP_Type h_opt = h_var * std::pow(TOL/local_error, 1./(order+1));

        if (h_opt < h_var)
          { // (3) Time step is rejected; repeat step with optimal value.
            h_var = h_opt;
            misfires++;
            continue;
          }
        else
          { // (4) Time step was accepted.
            t += h_var;
            u += h_var * inc_u;

            EOS_Method::save_step(t, u);
            step_sizes.emplace_back(h_var);

            // Set time step for next iteration.
            h_var = h_opt;

            // Avoid rounding errors (Remark 2.4.3)
            if (t + 1.1*h_var >= t_lim)
              h_var = t_lim - t;
          }
      }
  }

  void print_step_size(std::ostream &out)
  {
    for (size_t i = 0; i < step_sizes.size(); i++)
      {
        out << i << "\t" << step_sizes[i] << std::endl;
      }
  }

private:
  // EOS_Method::f;
  // EOS_Method::t0;
  // EOS_Method::u0;
  // EOS_Method::steps;
  // EOS_Method::timepoints;
  // EOS_Method::uapprox;

  // Butcher tableau
  dealii::LAPACKFullMatrix<FP_Type> A;
  dealii::Vector<FP_Type> b1, b2, c;

  // Adaptive step-size plot
  std::vector<FP_Type> step_sizes;

  // Markers
  bool embedded_method;
  size_t misfires;
};

#endif // RUNGE_KUTTA_H
