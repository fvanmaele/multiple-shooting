#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>

#include <deal.II/lac/full_matrix.h>

#include "../lac/vector_operators.h"
#include "tableau.h"
#include "eos_method.h"

class ERK : public EOS_Method
{
public:
  // In addition to the usual arguments, the Range_Kuuta method
  // requires a matrix A and two vectors b, c. The matrix A is stored
  // as a dense matrix, for simplicity and due to its small dimension.
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::FullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b1,
      dealii::Vector<FP_Type> _c) :
    EOS_Method(f, t0, u0), A(_A), b1(_b1), c(_c)
  {
    assert(c.size() == b1.size());
    embedded_method = false;
  }

  // Constructor for embedded methods
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::FullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b1,
      dealii::Vector<FP_Type> _b2,
      dealii::Vector<FP_Type> _c) :
    EOS_Method(f, t0, u0), A(_A), b1(_b1), b2(_b2), c(_c)
  {
    assert(c.size() == b1.size());
    assert(c.size() == b2.size());
    embedded_method = true;
  }

  dealii::Vector<FP_Type>
  k_increment(const FP_Type &t, const dealii::Vector<FP_Type> &u,
              const FP_Type &h, const dealii::Vector<FP_Type> &b)
  {
    size_t s = b.size();
    std::vector<dealii::Vector<FP_Type> > k(s);
    k.at(0) = f(t, u);

    for (size_t i = 1; i < s; i++)
      {
        dealii::Vector<FP_Type> g(u.size());

        for (size_t j = 0; j < i; j++)
          g += A(i,j) * k.at(j);

        k.at(i) = f(t + h*c[i], u + h*g);
        assert(k.at(i).size() == u.size());
      }

    dealii::Vector<FP_Type> result(u.size());

    for (size_t i = 0; i < s; i++)
      result += b(i) * k.at(i);
    return result;
  }

  virtual dealii::Vector<FP_Type>
  increment_function(const FP_Type &t, const dealii::Vector<FP_Type> &u,
                     const FP_Type &h) override
  {
    // Explicit method of higher order (constructor)
    return k_increment(t, u, h, b1);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  void iterate_with_ssc(const FP_Type t_lim, const FP_Type h0,
                        const FP_Type TOL, const size_t order)
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
        inc_u = u + h_var * k_increment(t, u, h_var, b1);
        inc_v = u + h_var * k_increment(t, u, h_var, b2);
        steps++;

        // (2) Compute optimal step size.
        // Use norm properties to save operations.
        FP_Type local_error = (inc_u - inc_v).l2_norm();
        FP_Type h_opt = h_var * std::pow(TOL/local_error, 1./(order+1));

        if (h_opt < h_var)
          { // (3) Time step is rejected; repeat step with optimal value.
            h_var = h_opt;
            misfires++;
            continue;
          }
        else
          { // (4) Time step was accepted.
            u  = inc_u;
            t += h_var;

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
        out << timepoints[i] << "\t" << step_sizes[i] << std::endl;
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
  dealii::FullMatrix<FP_Type> A;
  dealii::Vector<FP_Type> b1, b2, c;

  // Adaptive step-size plot
  std::vector<FP_Type> step_sizes;

  // Markers
  bool embedded_method;
  size_t misfires;
};

#endif // RUNGE_KUTTA_H
