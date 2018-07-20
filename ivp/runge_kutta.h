#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H
#include "../lac/matrix_operators.h"
#include "../lac/vector_operators.h"
#include "base_method.h"

class ERK : public IVP_Method
{
public:
  // In addition to the usual arguments, the Range_Kuuta method
  // requires a matrix A and two vectors b, c. The matrix A is stored
  // as a dense matrix, for simplicity and due to its small dimension.
  ERK(Functor &f, NumberType t0, dealii::Vector<NumberType> u0,
      dealii::FullMatrix<NumberType> _A,
      dealii::Vector<NumberType> _b,
      dealii::Vector<NumberType> _c,
      NumberType h = 1.0e-2) :
    IVP_Method(f, t0, u0, h), A(_A), b(_b), c(_c)
  {}

  // Constructor for embedded methods
  ERK(Functor &f, NumberType t0, dealii::Vector<NumberType> u0,
      dealii::FullMatrix<NumberType> _A,
      dealii::Vector<NumberType> _b,
      dealii::Vector<NumberType> _b2,
      dealii::Vector<NumberType> _c,
      NumberType h = 1.0e-2) :
    IVP_Method(f, t0, u0, h), A(_A), b(_b), b2(_b2), c(_c)
  {}

  // Overload for embedded methods
  void iteration_step(dealii::Vector<NumberType> &u1,
                      dealii::Vector<NumberType> &u2,
                      NumberType &t, const NumberType &h)
  {
    ;
  }

  virtual void iteration_step(dealii::Vector<NumberType> &u,
                              NumberType &t, const NumberType &h) override
  {
    // See lecture notes, p.29
    assert(c.size() == b.size());
    size_t s = c.size();

    // Store values k_i for evaluation of u
    std::vector<dealii::Vector<NumberType> > k;

    // Optional, see Remark 2.3.4
    std::vector<dealii::Vector<NumberType> > g;

    // Compute k_1 ... k_s
    for (size_t i = 0; i < s; i++)
      {
        // Compute sum(a_ij * k_j)
        dealii::Vector<NumberType> g_sum(u.size());

        // Compute k_1 ... k_{i-1}
        // If i = 0, this loop is not executed, i.e. g_sum remains 0.
        for (size_t j = 0; j < i; j++)
          {
            // While the upper half of the matrix A is already ignored,
            // entries in the lower half may still be zero.
            NumberType a_ij = matrix_element(A, i, j);
            if (a_ij > 0)
              g_sum += a_ij * k.at(j);
          }

        g.push_back(u + h*g_sum);
        k.push_back(f(t + h*c(i), g.at(i)));
      }

    // Compute new step
    dealii::Vector<NumberType> k_sum(u.size());
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
  dealii::FullMatrix<NumberType> A;
  dealii::Vector<NumberType> b, b2, c;
};

struct ERK_04
{
  const NumberType entries[16] = {
    0.0, 0.0, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0
  };
  const std::vector<NumberType> b = { 1./6, 2./6, 2./6, 1./6 };
  const std::vector<NumberType> c = { 0.0, 0.5, 0.5, 1.0 };
};

struct DOPRI
{
  const NumberType entries[49] = {
    0,            0,             0,            0,          0,             0,       0,
    1./5,         0,             0,            0,          0,             0,       0,
    3./40,        9./40,         0,            0,          0,             0,       0,
    44./45,       -56./15,       32./9,        0,          0,             0,       0,
    19372./6561,  -25360./2187,  64448./6561,  -212./729,  0,             0,       0,
    9017./3168,   -355./33,      46732./5247,  49./176,    -5103./18656,  0,       0,
    35./384,      0,             500./1113,    125./192,   -2187./6784,   11./84,  0
  };

  const std::vector<NumberType> b1 = { 35./384, 0, 500./1113, 125./192, -2187./6784, 11./84, 0 };
  const std::vector<NumberType> b2 = { 5179./57600, 0, 7571./16695, 393./640, -92097./339200, 187./2100, 1./40 };
  const std::vector<NumberType> c  = { 0, 1./5, 3./10, 4./5, 8./9, 1, 1 };
};

class ERK_Test_O4 : public IVP_Method
{
public:
  using IVP_Method::IVP_Method;

  virtual void iteration_step(dealii::Vector<NumberType> &u, NumberType &t, const NumberType &h) override
  {
    dealii::Vector<NumberType> k1 = f(t, u);
    dealii::Vector<NumberType> k2 = f(t + 0.5*h, u + 0.5*h*k1);
    dealii::Vector<NumberType> k3 = f(t + 0.5*h, u + 0.5*h*k2);
    dealii::Vector<NumberType> k4 = f(t + h, u + h*k3);

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
