#ifndef FORWARD_AD_H
#define FORWARD_AD_H
// Based on examples from:
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/ad_example.cpp

#include <cassert>
#include <vector>

#include <Sacado.hpp>

#include "types.h"
#include "../lac/lac_types.h"

// Forward AD class (dynamic allocation)
typedef Sacado::Fad::DFad<FP_Type> NumberAD;

// Save FAD input values in std::vector.
// https://github.com/dealii/dealii/issues/6940
typedef std::vector<NumberAD> VectorAD;

/*! \class FAD_Setup
 *  \brief Class representing the time-dependent function with \f$n\f$ components
 *  \f$f:I\times\mathbb{R}^m \rightarrow \mathbb{R}^n\f$
 *
 * Both evaluation and automatic differentation (using the Sacado
 * package from the Trilinos library) are supported.
 *
 * For general information, see SESS 2007, E. Phipps:
 *
 * https://software.sandia.gov/SESS/past_seminars/111307_Phipps.html
 */
template <typename Callable>
class FAD_Setup
{
public:
  /*! \fn FAD_Setup
   * \brief Constructor. Takes an initial value for the number of components.
   */
  FAD_Setup(Callable _f, size_t _m, size_t _n) :
    f(_f), u_dim(_m), f_dim(_n), f_init(0), FAD_t(0), FAD_u(_m), FAD_y(_n)
  {}

  /*! \fn FAD_Setup
   * \brief Constructor for equidimensional problems
   * \f$f:I\times\mathbb{R}^d \rightarrow \mathbb{R}^d\f$.
   */
  FAD_Setup(Callable _f, size_t _d) :
    f(_f), u_dim(_d), f_dim(_d), f_init(0), FAD_t(0), FAD_u(_d), FAD_y(_d)
  {}

  /*! \fn init
   *  \brief Evaluate function \f$(t, u)\f$ on AD variables.
   * Results may be retrieved using the \c value() and \c diff() methods.
   *
   * \f$u_1,\cdots,u_m\f$ are set as independent variables.
   */
  void init(FP_Type t, const VectorD2 &u)
  {
    if (u.size() != u_dim)
      throw std::domain_error("u is not in domain of f");

    // Passive variable
    FAD_t = t;

    // Analytic derivative with respect to u
    for (size_t i = 0; i < u_dim; i++)
      {
        FAD_u.at(i) = u(i);
        FAD_u.at(i).diff(i, u_dim);
      }

    FAD_y = f(FAD_t, FAD_u);
    f_init = true;

    if (FAD_y.size() != f_dim)
      throw std::range_error("y is not in range of f");
  }

  /*! \fn value
   *  \brief Return the value \f$(t, u)\f$ as a \c dealii vector.
   */
  VectorD2 value() const
  {
    if (!f_init)
      throw std::invalid_argument("FAD must be initialized");
    VectorD2 y(f_dim);

    for (size_t i = 0; i < FAD_y.size(); i++)
      y(i) = FAD_y.at(i).val();

    return y;
  }

  /*! \fn diff
   *  \brief Return the partial derivatives
   * \f$\frac{\partial f}{\partial u_1},\cdots,\frac{\partial f}{\partial u_m}\f$
   * as a \c dealii matrix.
   */
  MatrixD2 diff() const
  {
    if (!f_init)
      throw std::invalid_argument("FAD must be initialized");
    MatrixD2 J(f_dim, u_dim);

    for (size_t i = 0; i < f_dim; i++)
      if (FAD_y.at(i).hasFastAccess())
        for (size_t j = 0; j < u_dim; j++)
          J.set(i, j, FAD_y.at(i).fastAccessDx(j));
      else
        for (size_t j = 0; j < u_dim; j++)
          J.set(i, j, FAD_y.at(i).dx(j));

    return J;
  }

private:
  Callable f;
  size_t u_dim, f_dim;
  bool f_init;

  NumberAD FAD_t;
  VectorAD FAD_u;
  VectorAD FAD_y;
};

/*! \class FAD_tWrapper
 *  \brief Class for equidimensional, time-dependent problems
 * \f$f:I\times\mathbb{R}^d \rightarrow \mathbb{R}^d\f$ using AD.
 *
 * Adapter for \c TimeDivFunctor.
 */
template <typename Callable>
class FAD_tWrapper : public TimeDivFunctor
{
public:
  FAD_tWrapper(Callable f, size_t dim) :
    TimeDivFunctor(dim), F(f, dim)
  {}

  virtual VectorD2
  operator()(FP_Type t, const VectorD2 &u) override
  {
    F.init(t, u);
    return F.value();
  }

  virtual MatrixD2
  diff(FP_Type t, const VectorD2 &u) override
  {
    F.init(t, u);
    return F.diff();
  }

private:
  FAD_Setup<Callable> F;
  FP_Type t;
};

/*! \class FAD_cWrapper
 * \brief Class for equidimensional problems
 * \f$f:\mathbb{R}^d \rightarrow \mathbb{R}^d\f$ using AD.
 *
 * Adapter for \c DivFunctor.
 */
template <typename Callable>
class FAD_cWrapper : public DivFunctor
{
public:
  FAD_cWrapper(Callable f, size_t dim, FP_Type _t = 0) :
    DivFunctor(dim), F(f, dim), t(_t)
  {}

  virtual VectorD2
  operator()(const VectorD2 &u) override
  {
    F.init(t, u);
    return F.value();
  }

  virtual MatrixD2
  diff(const VectorD2 &u) override
  {
    F.init(t, u);
    return F.diff();
  }

private:
  FAD_Setup<Callable> F;
  FP_Type t;
};

#endif // FORWARD_AD_H
