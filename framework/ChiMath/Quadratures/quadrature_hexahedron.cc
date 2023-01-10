#include "quadrature_hexahedron.h"

#include "quadrature_gausslegendre.h"

//###################################################################
/**Initialzes a set of points for a quadrature integration over
 * the volume of a hexahedron.*/
chi_math::QuadratureHexahedron::QuadratureHexahedron(QuadratureOrder order) :
  Quadrature(order)
{
  QuadratureGaussLegendre legendre(order);

  size_t N = legendre.qpoints.size();

  qpoints.resize(N*N*N);
  weights.resize(N*N*N);

  unsigned int q=0;

  for (unsigned int k=0; k<N; k++)
    for (unsigned int j=0; j<N; j++)
      for (unsigned int i=0; i<N; i++)
      {
        qpoints[q](0) = legendre.qpoints[i](0);
        qpoints[q](1) = legendre.qpoints[j](0);
        qpoints[q](2) = legendre.qpoints[k](0);

        weights[q] = legendre.weights[i] *
                     legendre.weights[j] *
                     legendre.weights[k];

        q++;
      }
}