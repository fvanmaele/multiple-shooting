#ifndef LINEAR_H
#define LINEAR_H

#include "shooting.h"

// Class to represent a 2-dimensional separated linear BVP, with
// boundary conditions on the first component.
template <typename DiffMethod>
class SimpleBVP
{
public:
  SimpleBVP(TimeFunctor &_f, FP_Type _a, FP_Type _b, dealii::Vector<FP_Type> _c) :
    f(_f), a(_a), b(_b), c(_c)
  {
    std::array<FP_Type, 4> _A = { 1, 0,
                                  0, 0 };
    std::array<FP_Type, 4> _B = { 0, 0,
                                  1, 0 };
    // construct matrix
    A = dealii::FullMatrix<FP_Type>(2, 2, _A.data());
    B = dealii::FullMatrix<FP_Type>(2, 2, _B.data());
  }

  // Apply single shooting method applied to separated BVP.
  // As the Newton method may converge to different roots, take a vector
  // of starting values instead of a single entry.
  void single_shooting(const std::vector<dealii::Vector<FP_Type> > &start)
  {
    assert(start.size());
    DiffMethod F(f, a, b, A, B, c);
    Newton N(F, 2);

    for (auto &s : start)
      {
        if (std::is_same<DiffMethod, SF_External>::value)
          N.iterate(s);
        else
          N.iterate_broyden(s);
      }
  }

  void shooting_graph(size_t dim, const std::vector<dealii::Vector<FP_Type> > &range,
                      std::ofstream &output_file)
  {
    if (dim > 2)
      throw std::logic_error("Function not implemented");
    DiffMethod F(f, a, b, A, B, c);

    for (auto &s : range)
      {
        dealii::Vector<FP_Type> diff = F(s);
        switch (dim)
          {
          case 1 :
            output_file << s[0]
                << "\t" << diff[0] << std::endl;
            break;
          case 2 :
            output_file << s[0]
                << "\t" << s[1]
                << "\t" << diff[0]
                << "\t" << diff[1] << std::endl;
            break;
          }
      }
  }

private:
  TimeFunctor &f;
  FP_Type a, b;
  dealii::Vector<FP_Type> c;
  dealii::FullMatrix<FP_Type> A, B;
};

#endif // LINEAR_H
