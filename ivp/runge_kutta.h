#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>
#include <utility>

#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"
#include "eos_method.h"
#include "tableau.h"

template <typename ButcherTableau>
class ERK : public EOS_Method
{
public:
  // Constructor for explicit and embedded methods
  ERK(TimeFunctor &f, FP_Type t0, dealii::Vector<FP_Type> u0) :
    EOS_Method(f, t0, u0)
  {
    ButcherTableau BTab;

    A  = dealii::FullMatrix(BTab.n, BTab.n, BTab.A.data());
    b1 = dealii::Vector<FP_Type>(BTab.b_high.begin(), BTab.b_high.end());
    c  = dealii::Vector<FP_Type>(BTab.c.begin(), BTab.c.end());
    order = BTab.p;

    if (BTab.b_low.size())
      {
        embedded_method = true;
        b2 = dealii::Vector<FP_Type>(BTab.b_low.begin(), BTab.b_low.end());
        assert(b1.size() == b2.size());
      }
    else
      {
        embedded_method = false;
        b2 = dealii::Vector<FP_Type>();
      }
  }

  dealii::Vector<FP_Type>
  k_increment(FP_Type t, const dealii::Vector<FP_Type> &u,
              FP_Type h, const dealii::Vector<FP_Type> &b)
  {
    size_t s = b.size();
    size_t n = u.size();
    std::vector<dealii::Vector<FP_Type> > k(s);
    k.at(0) = f(t, u);

    for (size_t i = 1; i < s; i++)
      {
        dealii::Vector<FP_Type> g(u.size());

        for (size_t j = 0; j < i; j++)
          g += A(i, j) * k.at(j);

        k.at(i) = f(t + h*c[i], u + h*g);
      }

    dealii::Vector<FP_Type> u_inc(n);

    for (size_t i = 0; i < s; i++)
      u_inc += b(i) * k.at(i);

    return u_inc;
  }


  std::pair<dealii::Vector<FP_Type>, dealii::FullMatrix<FP_Type> >
  k_variational(FP_Type t, const dealii::Vector<FP_Type> &u,
                FP_Type h, const dealii::Vector<FP_Type> &b,
                const dealii::FullMatrix<FP_Type> &Y,
                TimeDivFunctor *F)
  {
    size_t s = b.size();
    size_t n = u.size();
    assert(Y.m() == n);
    assert(Y.n() == n);

    // Increments k_1, ..., k_s for solution u(t)
    std::vector<dealii::Vector<FP_Type> > k(s);
    k.at(0) = (*F)(t, u);

    // Increments K_1, ..., K_s for solution Y(t)
    std::vector<dealii::FullMatrix<FP_Type> > K(s);
    K.at(0) = (*F).diff(t, u) * Y;

    for (size_t i = 1; i < s; i++)
      {
        dealii::Vector<FP_Type> g(n);
        dealii::FullMatrix<FP_Type> G(n);

        for (size_t j = 0; j < i; j++)
          {
            g += A(i, j) * k.at(j);
            G.add(1, A(i, j) * K.at(j));
          }

        k.at(i) = (*F)(t + h*c[i], u + h*g);
        K.at(i) = (*F).diff(t + h*c[i], h*g) * (Y + h*G);
      }

    dealii::Vector<FP_Type> inc_u(n);
    dealii::FullMatrix<FP_Type> inc_Y(n, n);

    for (size_t i = 0; i < s; i++)
      {
        inc_u += b(i) * k.at(i);
        inc_Y.add(1, b(i) * K.at(i));
      }
    return std::make_pair(inc_u, inc_Y);
  }

  virtual dealii::Vector<FP_Type>
  increment_function(FP_Type t, const dealii::Vector<FP_Type> &u,
                     FP_Type h) override
  {
    // Explicit method of higher order
    return k_increment(t, u, h, b1);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  dealii::FullMatrix<FP_Type>
  fund_matrix() const
  {
    return Y;
  }

  void iterate_with_ssc(const FP_Type t_lim, const FP_Type h0,
                        const FP_Type TOL, bool fundamental_matrix)
  {
    assert(embedded_method);
    FP_Type t = t0;     // start time
    FP_Type h_var = h0; // initial step size
    size_t n = u0.size();

    TimeDivFunctor* f_diff = dynamic_cast<TimeDivFunctor*>(&f);
    if (fundamental_matrix && f_diff == nullptr)
      throw std::invalid_argument("right-hand side is not differentiable");

    // Dynamic allocation, declare outside loop
    dealii::Vector<FP_Type> u = u0;
    dealii::Vector<FP_Type> inc_u(n);
    dealii::Vector<FP_Type> inc_v(n);

    // Init output variables
    EOS_Method::reset();
    steps = 0;
    misfires = 0;

    // Initial value for fundamental matrix Y(t; t0)
    Y = dealii::IdentityMatrix(n);
    dealii::FullMatrix<FP_Type> inc_Y(n, n);

    // Algorithm 2.4.2
    // The step-size is controlled only by the IVP, not the variational equation,
    // to avoid recomputing the Jacobian on a rejected time step.
    while (t_lim - t > 0)
      { // (1) Candidates for u_k, v_k with time step h_k
        if (fundamental_matrix)
          {
            auto U = k_variational(t, u, h_var, b1, Y, f_diff);

            inc_u = u + h_var * U.first;
            inc_v = u + h_var * k_increment(t, u, h_var, b2);
            inc_Y = Y + h_var * U.second;
          }
        else
          {
            inc_u = u + h_var * k_increment(t, u, h_var, b1);
            inc_v = u + h_var * k_increment(t, u, h_var, b2);
          }
        steps++;

        // (2) Compute optimal step size.
        FP_Type local_error = (inc_u - inc_v).l2_norm();
        FP_Type h_opt = h_var * std::pow(TOL / local_error, 1. / order);

        // (3) Time step is rejected; repeat step with optimal value.
        if (h_opt < h_var)
          {
            h_var = h_opt;
            misfires++;
            continue;
          }
        else
          { // (4) Time step was accepted.
            u = inc_u;
            Y = inc_Y;
            t += h_var;

            // Save intermediary solutions to IVP. (In the variational equation,
            // only the last step Y_n is stored.)
            EOS_Method::save_step(t, u);
            step_sizes.push_back(h_var);

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
  size_t order;

  // Variational equation
  dealii::FullMatrix<FP_Type> Y;

  // Adaptive step-size plot
  std::vector<FP_Type> step_sizes;

  // Markers
  bool embedded_method;
  size_t misfires;
};

#endif // RUNGE_KUTTA_H
