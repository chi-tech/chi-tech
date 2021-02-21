#ifndef QUADRATURE_GAUSS_CHEBYSHEV_H
#define QUADRATURE_GAUSS_CHEBYSHEV_H

#include "quadrature.h"

namespace chi_math
{
  class QuadratureGaussChebyshev;
}

//######################################################### Class Def
/**Gauss-Chebyshev quadrature.*/
class chi_math::QuadratureGaussChebyshev : public chi_math::Quadrature
{
public:
  QuadratureGaussChebyshev() :
    chi_math::Quadrature(QuadratureOrder::CONSTANT)
  {}
  //01
  void Initialize(int N, bool verbose=false);

};

#endif