#include "quadrature_conical.h"

#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/Jacobi/quadrature_jacobi.h"

//###################################################################
/**Initialize conical quadrature for a triangle.*/
void chi_math::QuadratureConical::Initialize_Conical_Product_Tri()
{
  QuadratureGaussLegendre legendre((QuadratureOrder)std::floor(((int)order+1)/2.0));
  QuadratureJacobi        jacobiA(order, /*alpha=*/1, /*beta=*/0);

  legendre.Scale({-1.0, 1.0},
                 { 0.0, 1.0});

  const size_t np = legendre.qpoints.size();

  qpoints.resize(np * np);
  weights.resize(np * np);

  // Compute the conical product
  unsigned int gp = 0;
  for (unsigned int i=0; i<np; i++)
    for (unsigned int j=0; j<np; j++)
    {
      qpoints[gp](0) = jacobiA.qpoints[j](0);           //s[j];
      qpoints[gp](1) = legendre.qpoints[i](0) *
                       (1.-jacobiA.qpoints[j](0));      //r[i]*(1.-s[j]);
      weights[gp]    = legendre.weights[i] *
                       jacobiA.weights[j];              //A[i]*B[j];
      gp++;
    }
}

//###################################################################
/**Initialize conical quadrature for a tetrahedron.*/
void chi_math::QuadratureConical::Initialize_Conical_Product_Tet()
{
  QuadratureGaussLegendre legendre((QuadratureOrder)std::floor(((int)order+1)/2.0));
  QuadratureJacobi        jacobiA(order, /*alpha=*/1, /*beta=*/0);
  QuadratureJacobi        jacobiB(order, /*alpha=*/2, /*beta=*/0);

  legendre.Scale({-1.0, 1.0},
                 { 0.0, 1.0});

  const size_t np = legendre.qpoints.size();

  qpoints.resize(np * np * np);
  weights.resize(np*np*np);

  double weight_sum = 0.0;
  unsigned int gp = 0;
  for (unsigned int i=0; i<np; i++)
    for (unsigned int j=0; j<np; j++)
      for (unsigned int k=0; k<np; k++)
      {
        qpoints[gp](0) = jacobiB.qpoints[k](0);       //t[k];
        qpoints[gp](1) = jacobiA.qpoints [j](0) *
                         (1.-jacobiB.qpoints[k](0));  //s[j]*(1.-t[k]);
        qpoints[gp](2) = legendre.qpoints[i](0) *
                         (1.-jacobiA.qpoints[j](0)) *
                         (1.-jacobiB.qpoints[k](0));  //r[i]*(1.-s[j])*(1.-t[k]);
        weights[gp]   = legendre.weights[i]*
                        jacobiA.weights [j]*
                        jacobiB.weights [k];              //A[i]*B[j]*C[k];
        weight_sum += weights[gp];
        gp++;
      }

  double w_scale = (1.0/6.0)/weight_sum;
  for (auto& v : weights)
    v *= w_scale;
}