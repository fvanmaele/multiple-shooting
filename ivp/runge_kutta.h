#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>
#include <utility>

#include "../lac/vector_operators.h"
#include "../lac/matrix_operators.h"
#include "eos_method.h"
#include "tableau.h"

template <typename ButcherTableau>
class ERK : public OneStepMethod
{
public:
  // Constructor for explicit and embedded methods
  ERK(TimeFunctor &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      Curve *u = nullptr) :
    OneStepMethod(f, t0, u0, u)
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
  k_increment(FP_Type t, const dealii::Vector<FP_Type> &y,
              FP_Type h, const dealii::Vector<FP_Type> &b)
  {
    size_t s = b.size();
    size_t n = y.size();
    std::vector<dealii::Vector<FP_Type> > k(s);
    k.at(0) = f(t, y);

    for (size_t i = 1; i < s; i++)
      {
        dealii::Vector<FP_Type> g(y.size());

        for (size_t j = 0; j < i; j++)
          g += A(i, j) * k.at(j);

        k.at(i) = f(t + h*c[i], y + h*g);
      }

    dealii::Vector<FP_Type> u_inc(n);

    for (size_t i = 0; i < s; i++)
      u_inc += b(i) * k.at(i);

    return u_inc;
  }


  std::pair<dealii::Vector<FP_Type>, dealii::FullMatrix<FP_Type> >
  k_variational(FP_Type t, const dealii::Vector<FP_Type> &y,
                FP_Type h, const dealii::Vector<FP_Type> &b,
                const dealii::FullMatrix<FP_Type> &Y,
                TimeDivFunctor *F)
  {
    size_t s = b.size();
    size_t n = y.size();
    assert(Y.m() == n);
    assert(Y.n() == n);

    // Increments k_1, ..., k_s for solution y(t)
    std::vector<dealii::Vector<FP_Type> > k(s);
    k.at(0) = (*F)(t, y);

    // Increments K_1, ..., K_s for solution Y(t)
    std::vector<dealii::FullMatrix<FP_Type> > K(s);
    K.at(0) = (*F).diff(t, y) * Y;

    for (size_t i = 1; i < s; i++)
      {
        dealii::Vector<FP_Type> g(n);
        dealii::FullMatrix<FP_Type> G(n);

        for (size_t j = 0; j < i; j++)
          {
            g += A(i, j) * k.at(j);
            G.add(1, A(i, j) * K.at(j));
          }

        k.at(i) = (*F)(t + h*c[i], y + h*g);
        K.at(i) = (*F).diff(t + h*c[i], y + h*g) * (Y + h*G);
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
  increment_function(FP_Type t, const dealii::Vector<FP_Type> &y,
                     FP_Type h) override
  {
    // Explicit method of higher order
    return k_increment(t, y, h, b1);
  }

  size_t n_misfires() const
  {
    return misfires;
  }

  dealii::FullMatrix<FP_Type>
  fund_matrix() const
  {
    return Yn;
  }

  void iterate_with_ssc(FP_Type t_lim, FP_Type h0, FP_Type TOL,
                        bool fundamental_matrix = false,
                        FP_Type C = 2)
  {
    assert(embedded_method);
    FP_Type t = t0;     // start time
    FP_Type h_var = h0; // initial step size
    size_t n = u0.size();

    TimeDivFunctor* f_diff = dynamic_cast<TimeDivFunctor*>(&f);
    if (fundamental_matrix && f_diff == nullptr)
      throw std::invalid_argument("right-hand side is not differentiable");

    // Dynamic allocation, declare outside loop
    dealii::Vector<FP_Type> y = u0;
    dealii::Vector<FP_Type> inc_y1(n);
    dealii::Vector<FP_Type> inc_y2(n);

    // Init output variables
    OneStepMethod::reset();
    steps = 0;
    misfires = 0;

    // Initial value for fundamental matrix Yn(t; t0)
    Yn = dealii::IdentityMatrix(n);
    dealii::FullMatrix<FP_Type> inc_Yn(n, n);

    // Algorithm 2.4.2
    // The step-size is controlled only by the IVP, not the variational equation,
    // to avoid recomputing the Jacobian on a rejected time step.
    while (t_lim - t > 0)
      { // (1) Candidates for u_k, v_k with time step h_k
        if (fundamental_matrix)
          {
            auto U = k_variational(t, y, h_var, b1, Yn, f_diff);

            inc_y1 = y + h_var * U.first;
            inc_y2 = y + h_var * k_increment(t, y, h_var, b2);
            inc_Yn = Yn + h_var * U.second;
          }
        else
          {
            inc_y1 = y + h_var * k_increment(t, y, h_var, b1);
            inc_y2 = y + h_var * k_increment(t, y, h_var, b2);
          }
        steps++;

        if (sol_is_nan(inc_y1))
          throw std::overflow_error("global error too large (NaN)");

        // (2) Compute optimal step size.
        FP_Type local_error = (inc_y1 - inc_y2).l2_norm();
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
            y  = inc_y1;
            Yn = inc_Yn;
            t += h_var;

            // Save intermediary solutions to IVP. (In the variational equation,
            // only the last step Y_n is stored.)
            OneStepMethod::save_step(t, y);
            step_sizes.push_back(h_var);

            // Set time step for next iteration.
            h_var = h_opt;

            // (*) Guarantee to hit right interval end (Rem. 2.4.3)
            if (t + 1.1*h_var >= t_lim)
              h_var = t_lim - t;

            // Comparison to exact solution
            if (u != nullptr)
              if (y.l2_norm() >= C*(*u)(t).l2_norm())
                throw std::domain_error("global error too large");
          }
      }

    if (timepoints.back() != t_lim)
      std::cerr << "warning: time step outside interval boundary ("
                << timepoints.back() << "; [" << t0 << ", " << t_lim << "])"
                << std::endl;
  }

  void print_step_size(std::ostream &out)
  {
    for (size_t i = 0; i < step_sizes.size(); i++)
      {
        out << timepoints[i] << "\t" << step_sizes[i] << std::endl;
      }
  }

private:
  // OneStepMethod::f;
  // OneStepMethod::u;
  // OneStepMethod::t0;
  // OneStepMethod::u0;
  // OneStepMethod::steps;
  // OneStepMethod::timepoints;
  // OneStepMethod::uapprox;

  // Butcher tableau
  dealii::FullMatrix<FP_Type> A;
  dealii::Vector<FP_Type> b1, b2, c;
  size_t order;

  // Variational equation
  dealii::FullMatrix<FP_Type> Yn;

  // Adaptive step-size plot
  std::vector<FP_Type> step_sizes;

  // Markers
  bool embedded_method;
  size_t misfires;
};

#endif // RUNGE_KUTTA_H
