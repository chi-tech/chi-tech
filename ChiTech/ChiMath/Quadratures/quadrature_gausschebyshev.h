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
  explicit
  QuadratureGaussChebyshev(QuadratureOrder in_order,
                           bool verbose=false);

  explicit
  QuadratureGaussChebyshev(unsigned int N,
                           bool verbose=false);

private:
  void Initialize(unsigned int N, bool verbose);
};

#endif