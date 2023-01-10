#include "quadrature_quadrilateral.h"

#include "quadrature_gausslegendre.h"

//###################################################################
/**Initialzes a set of points for a quadrature integration over
 * the volume of a quadrilateral.*/
chi_math::QuadratureQuadrilateral::QuadratureQuadrilateral(QuadratureOrder order) :
  Quadrature(order)
{
  QuadratureGaussLegendre legendre(order);

  size_t N = legendre.qpoints.size();

  qpoints.resize(N*N);
  weights.resize(N*N);

  unsigned int q=0;

  for (unsigned int j=0; j<N; j++)
    for (unsigned int i=0; i<N; i++)
    {
      qpoints[q](0) = legendre.qpoints[i](0);
      qpoints[q](1) = legendre.qpoints[j](0);

      weights[q] = legendre.weights[i] *
                   legendre.weights[j];

      q++;
    }
}