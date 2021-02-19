#ifndef QUADRATURE_JACOBI_H
#define QUADRATURE_JACOBI_H

#include "ChiMath/Quadratures/quadrature.h"

namespace chi_math
{


//###################################################################
/**Jacobi quadrature.*/
class QuadratureJacobi : public chi_math::Quadrature
{
private:
  const unsigned int m_alpha;
  const unsigned int m_beta;
public:
  QuadratureJacobi(QuadratureOrder order,
                   unsigned int a,
                   unsigned int b) :
                   m_alpha(a),
                   m_beta(b)
  {
    Initialize(order);
  }

private:
  void Initialize(QuadratureOrder order);
};



}//namespace chi_math

#endif //QUADRATURE_JACOBI_H