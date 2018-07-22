#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "../lac/vector_operators.h"
#include "base_method.h"
#include "butcher_tableau.h"

class ERK : public IVP_Method
{
public:
  // In addition to the usual arguments, the Range_Kuuta method
  // requires a matrix A and two vectors b, c. The matrix A is stored
  // as a dense matrix, for simplicity and due to its small dimension.
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::LAPACKFullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b,
      dealii::Vector<FP_Type> _c,
      FP_Type h = 1.0e-2) :
    IVP_Method(f, t0, u0, h), A(_A), b(_b), c(_c)
  {}

  // Constructor for embedded methods
  ERK(RHS &f, FP_Type t0, dealii::Vector<FP_Type> u0,
      dealii::LAPACKFullMatrix<FP_Type> _A,
      dealii::Vector<FP_Type> _b,
      dealii::Vector<FP_Type> _b2,
      dealii::Vector<FP_Type> _c,
      FP_Type h = 1.0e-2) :
    IVP_Method(f, t0, u0, h), A(_A), b(_b), b2(_b2), c(_c)
  {}

  // Overload for embedded methods
  void iteration_step(dealii::Vector<FP_Type> &u1,
                      dealii::Vector<FP_Type> &u2,
                      FP_Type &t, const FP_Type &h)
  {
    ;
  }

  virtual void
  iteration_step(dealii::Vector<FP_Type> &u,
                 FP_Type &t, const FP_Type &h) override
  {
    // See lecture notes, p.29
    assert(c.size() == b.size());
    size_t s = c.size();

    // Store values k_i for evaluation of u
    std::vector<dealii::Vector<FP_Type> > k;

    // Optional, see Remark 2.3.4
    std::vector<dealii::Vector<FP_Type> > g;

    // Compute k_1 ... k_s
    for (size_t i = 0; i < s; i++)
      {
        // Compute sum(a_ij * k_j)
        dealii::Vector<FP_Type> g_sum(u.size());

        // Compute k_1 ... k_{i-1}
        // If i = 0, this loop is not executed, i.e. g_sum remains 0.
        for (size_t j = 0; j < i; j++)
          {
            // While the upper half of the matrix A is already ignored,
            // entries in the lower half may still be zero.
            FP_Type a_ij = A(i, j);
            if (a_ij > 0)
              g_sum += a_ij * k.at(j);
          }

        g.push_back(u + h*g_sum);
        k.push_back(f(t + h*c(i), g.at(i)));
      }

    // Compute new step
    dealii::Vector<FP_Type> k_sum(u.size());
    for (size_t i = 0; i < s; i++)
      {
        k_sum += b(i) * k.at(i);
      }
    u += h * k_sum;
    t += h;
  }

private:
  // IVP_Method::f;
  // IVP_Method::t0;
  // IVP_Method::u0;
  // IVP_Method::h;
  // IVP_Method::nsteps;

  dealii::LAPACKFullMatrix<FP_Type> A;
  dealii::Vector<FP_Type> b, b2, c;
};

class ERK_Test_O4 : public IVP_Method
{
public:
  using IVP_Method::IVP_Method;

  virtual void iteration_step(dealii::Vector<FP_Type> &u, FP_Type &t, const FP_Type &h) override
  {
    dealii::Vector<FP_Type> k1 = f(t, u);
    dealii::Vector<FP_Type> k2 = f(t + 0.5*h, u + 0.5*h*k1);
    dealii::Vector<FP_Type> k3 = f(t + 0.5*h, u + 0.5*h*k2);
    dealii::Vector<FP_Type> k4 = f(t + h, u + h*k3);

    u += h * (1./6*k1 + 2./6*k2 + 2./6*k3 + 1./6*k4);
    t += h;
  }

private:
  // IVP_Method::f;
  // IVP_Method::t0;
  // IVP_Method::u0;
  // IVP_Method::h;
  // IVP_Method::nsteps;
};

#endif // RUNGE_KUTTA_H
