#ifndef TABLEAU_H
#define TABLEAU_H

#include <array>

#include "../base/types.h"

// Classic Range-Kutta method
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

// Cash-Karp method
struct KARP
{
  const size_t n = 6;
  const size_t p = 5;

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

// Dormand-Prince method
struct DOPRI
{
  const size_t n = 7;
  const size_t p = 5;

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

struct DOPRI87
{
  const size_t n = 13;
  const size_t p = 8;

  const std::array<FP_Type, 13> c =  {
    0., 1./18, 1./12, 1./8, 5./16, 3./8, 59./400, 93./200, 5490023248./9719169821, 13./20, 1201146811./1299019798, 1., 1.
  };

  const std::array<FP_Type, 169> A = {
    0,                       0,      0,        0,                          0,                       0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    1./18,                   0,      0,        0,                          0,                       0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    1./48,                   1./16,  0,        0,                          0,                       0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    1./32,                   0,      3./32,    0,                          0,                       0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    5./16,                   0,      -75./64,  75./64,                     0,                       0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    3./80,                   0,      0,        3./16,                      3./20,                   0,                         0,                         0,                         0,                         0,                       0,                      0,  0,
    29443841./614563906,     0,      0,        77736538./692538347,        -28693883./1125000000,   23124283./1800000000,      0,                         0,                         0,                         0,                       0,                      0,  0,
    16016141./946692911,     0,      0,        61564180./158732637,        22789713./633445777,     545815736./2771057229,     -180193667./1043307555,    0,                         0,                         0,                       0,                      0,  0,
    39632708./573591083,     0,      0,        -433636366./683701615,      -421739975./2616292301,  100302831./723423059,      790204164./839813087,      800635310./3783071287,     0,                         0,                       0,                      0,  0,
    246121993./1340847787,   0,      0,        -37695042795./15268766246,  -309121744./1061227803,  -12992083./490766935,      6005943493./2108947869,    393006217./1396673457,     123872331./1001029789,     0,                       0,                      0,  0,
    -1028468189./846180014,  0,      0,        8478235783./508512852,      1311729495./1432422823,  -10304129995./1701304382,  -48777925059./3047939560,  15336726248./1032824649,   -45442868181./3398467696,  3065993473./597172653,   0,                      0,  0,
    185892177./718116043,    0,      0,        -3185094517./667107341,     -477755414./1098053517,  -703635378./230739211,     5731566787./1027545527,    5232866602./850066563,     -4093664535./808688257,    3962137247./1805957418,  65686358./487910083,    0,  0,
    403863854./491063109,    0,      0,        -5068492393./434740067,     -411421997./543043805,   652783627./914296604,      11173962825./925320556,    -13158990841./6184727034,  3936647629./1978049680,    -160528059./685178525,   248638103./1413531060,  0,  0
  };

   const std::array<FP_Type, 13> b_high = {
     14005451./335480064, 0, 0, 0, 0, -59238493./1068277825, 181606767./758867731,   561292985./797845732,   -1041891430./1371343529,  760417239./1151165299, 118820643./751138087, -528747749./2220607170,  1./4
   };

   const std::array<FP_Type, 13> b_low = {
     13451932./455176623, 0, 0, 0, 0, -808719846./976000145, 1757004468./5645159321, 656045339./265891186,   -3867574721./1518517206,   465885868./322736535,  53011238./667516719, 2./45, 0
   };
};

#endif // TABLEAU_H
