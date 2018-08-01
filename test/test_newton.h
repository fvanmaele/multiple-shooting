#ifndef TEST_NEWTON_H
#define TEST_NEWTON_H
#include <cmath>
#include <vector>

#include "../algo/newton.h"
#include "../base/forward_ad.h"

struct FAD_Aut : public FAD_TdVecField
{
  virtual std::vector<FAD_Number>
  operator()(FAD_Number, const std::vector<FAD_Number> &x) override
  {
    std::vector<FAD_Number> y = {
      std::pow(x[0], 2) + std::pow(x[1], 2) - 1,
      std::pow(x[0], 2) - x[1]
    };
    return y;
  }
};

template <size_t dim>
class Functor_AD : public DivFunctor
{
public:
  Functor_AD(FAD_TdVecField &f) :
    F(f) // once differentiable
  {}

  virtual dealii::Vector<FP_Type>
  value(const dealii::Vector<FP_Type> &u) override
  {
    F.init(FP_Type(), u);
    return F.value();
  }

  virtual dealii::FullMatrix<FP_Type>
  jacobian(const dealii::Vector<FP_Type> &u) override
  {
    F.init(FP_Type(), u);
    return F.nabla_u();
  }

private:
  FAD_Wrapper<dim> F;
};

// Ortega, Scientific Computing, Tab. 4.3.2
void Test_Newton()
{
  FAD_Aut f;
  Functor_AD<2> F(f);

  // Starting value
  dealii::Vector<FP_Type> s(2);
  s[0] = 0.5;
  s[1] = 0.5;

  Newton_iterate(F, s);
}

#endif // TEST_NEWTON_H
