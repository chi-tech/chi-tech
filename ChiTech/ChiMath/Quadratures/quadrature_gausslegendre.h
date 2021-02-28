#ifndef QUADRATURE_GAUSS_LEGENDRE_H
#define QUADRATURE_GAUSS_LEGENDRE_H

#include "quadrature.h"
#include <stdio.h>

namespace chi_math
{
  class QuadratureGaussLegendre;
}

//######################################################### Class Def
/**Gauss-Legendre quadrature.*/
class chi_math::QuadratureGaussLegendre : public chi_math::Quadrature
{
public:
  //01
  explicit
  QuadratureGaussLegendre(QuadratureOrder in_order,
                          int maxiters=1000,
                          double tol=1.0e-12,
                          bool verbose=false);

  explicit
  QuadratureGaussLegendre(unsigned int N,
                          int maxiters=1000,
                          double tol=1.0e-12,
                          bool verbose=false);


  static std::vector<double> FindRoots(int N,
                                       int max_iters=1000,
                                       double tol=1.0e-12);
};

#endif