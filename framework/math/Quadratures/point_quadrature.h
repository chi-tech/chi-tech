#ifndef CHITECH_POINT_QUADRATURE_H
#define CHITECH_POINT_QUADRATURE_H

#include "quadrature.h"

namespace chi_math
{

/**Quadrate for a single point. Helps generalize quadrature based integration
* on 1D cell faces.*/
class PointQuadrature : public Quadrature
{
public:
  /**Quadrature order is not used and will always default to constant.*/
  explicit PointQuadrature(QuadratureOrder order = QuadratureOrder::CONSTANT);
};

}

#endif // CHITECH_POINT_QUADRATURE_H
