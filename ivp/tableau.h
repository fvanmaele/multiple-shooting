#ifndef TABLEAU_H
#define TABLEAU_H

#include <array>

#include "../base/types.h"

/*! \struct ERK_04
 * \brief Butcher Tableau for the classic Runge-Kutta method.
 */
struct ERK_04
{
  const size_t n = 4; // dimension
  const size_t p = 4; // order

  const std::array<FP_Type, 4> c = {
    0, 0.5, 0.5, 1.0
  };

  const std::array<FP_Type, 16> A = {
    0,    0,    0,    0,
    0.5,  0,    0,    0,
    0,    0.5,  0,    0,
    0,    0,    1.0,  0
  };

  const std::array<FP_Type, 4> b_high = {
    1./6, 2./6, 2./6, 1./6
  };

  const std::array<FP_Type, 0> b_low = {};
};

/*! \struct KARP
 * \brief Butcher Tableau for the Cash-Karp method of order 5(4).
 *
 * Ddetails on this method and its lower order variants are available
 * in the original paper:
 *
 * http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf
 */
struct KARP
{
  const size_t n = 6;
  const size_t p = 4; // low order

  const std::array<FP_Type, 6> c  = {
    0, 1./5, 3./10, 3./5, 1, 7./8
  };

  const std::array<FP_Type, 36> A = {
    0,            0,         0,           0,              0,          0,
    1./5,         0,         0,           0,              0,          0,
    3./40,        9./40,     0,           0,              0,          0,
    3./10,        -9./10,    6./5,        0,              0,          0,
    -11./54,      5./2,      -70./27,     35./27,         0,          0,
    1631./55296,  175./512,  575./13824,  44275./110592,  253./4096,  0
  };

  const std::array<FP_Type, 6> b_high = {
    37./378, 0, 250./621, 125./594, 0, 512./1771
  };

  const std::array<FP_Type, 6> b_low = {
    2825./27648, 0, 18575./48384, 13525./55296, 277./14336, 1./4
  };
};

/*! struct DOPRI54
 * \brief Butcher Tableau for the Dormand-Prince method of order 5(4).
 */
struct DOPRI54
{
  const size_t n = 7;
  const size_t p = 4; // low order

  const std::array<FP_Type, 7> c  = {
    0, 1./5, 3./10, 4./5, 8./9, 1., 1.
  };

  const std::array<FP_Type, 49> A = {
    0,            0,             0,            0,          0,             0,       0,
    1./5,         0,             0,            0,          0,             0,       0,
    3./40,        9./40,         0,            0,          0,             0,       0,
    44./45,       -56./15,       32./9,        0,          0,             0,       0,
    19372./6561,  -25360./2187,  64448./6561,  -212./729,  0,             0,       0,
    9017./3168,   -355./33,      46732./5247,  49./176,    -5103./18656,  0,       0,
    35./384,      0,             500./1113,    125./192,   -2187./6784,   11./84,  0
  };

  const std::array<FP_Type, 7> b_high = {
    35./384, 0, 500./1113, 125./192, -2187./6784, 11./84, 0
  };

  const std::array<FP_Type, 7> b_low = {
    5179./57600, 0, 7571./16695, 393./640, -92097./339200, 187./2100, 1./40
  };
};

#endif // TABLEAU_H
