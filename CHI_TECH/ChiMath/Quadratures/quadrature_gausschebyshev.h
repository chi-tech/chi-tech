#ifndef _quadrature_gauss_chebyshev_h
#define _quadrature_gauss_chebyshev_h

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
  //01
  void Initialize(int N, bool verbose=false);

};

#endif