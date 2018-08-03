#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>

#include <deal.II/lac/full_matrix.h>

#include "../lac/vector_operators.h"
#include "eos_method.h"
#include "tableau.h"

template <typename ButcherTableau>
class ERK : public EOS_Method
{
public:
  // Constructor for explicit and embedded methods
  ERK(tVecField &f, FP_Type t0, dealii::Vector<FP_Type> u0) :
    EOS_Method(f, t0, u0)
  {
    ButcherTableau BTab;

    // move assignment?
    A  = dealii::FullMatrix(BTab.n, BTab.n, BTab.A.data());
    b1 = dealii::Vector<FP_Type>(BTab.b_high.begin(), BTab.b_high.end());
    c  = dealii::Vector<FP_Type>(BTab.c.begin(), BTab.c.end());

    if (BTab.b_low.size())
      {
        embedded_method = true;
        b2 = dealii::Vector<FP_Type>(BTab.b_low.begin(), BTab.b_low.end());
      }
    else
      {
        embedded_method = false;
        b2 = dealii::Vector<FP_Type>();
      }
  }

  dealii::Vector<FP_Type>
  k_increment(FP_Type t, const dealii::Vector<FP_Type> &u,
              FP_Type h, const dealii::Vector<FP_Type> &b, tVecField &f)
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
  increment_function(FP_Type t, const dealii::Vector<FP_Type> &u,
                     FP_Type h, tVecField &f) override
  {
    // Explicit method of higher order (constructor)
    return k_increment(t, u, h, b1, f);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  void iterate_with_ssc(const FP_Type t_lim, const FP_Type h0,
                        const FP_Type TOL, const size_t order,
                        bool fundamental_matrix = false)
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
        inc_u = u + h_var * k_increment(t, u, h_var, b1, f);
        inc_v = u + h_var * k_increment(t, u, h_var, b2, f);
        steps++;

        // (2) Compute optimal step size.
        // Use norm properties to save operations.
        FP_Type local_error = (inc_u - inc_v).l2_norm();
        FP_Type h_opt = h_var * std::pow(TOL/local_error, 1./(order+1));

        // (3) Time step is rejected; repeat step with optimal value.
        if (h_opt < h_var)
          {
            h_var = h_opt;
            misfires++;
            continue;
          }
        // (4) Time step was accepted.
        else
          { // Take last values (t, u) to compute step in variational equation,
            // including the corresponding step size h_var.
            if (fundamental_matrix)
              // Note: the new step is written in-place, unlike the IVP solution.
              Y.add(h_var, fund_matrix_increment(t, u, h_var, Y, f));

            // Write new accepted values.
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
  // EOS_Method::Y;

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
