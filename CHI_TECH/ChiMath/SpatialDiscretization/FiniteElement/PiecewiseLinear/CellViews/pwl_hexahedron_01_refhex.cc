#include "pwl_hexahedron.h"

double HexahedronPWLFEValues::HexShape(unsigned int index,
                                       double xi,
                                       double eta,
                                       double zeta)
{
  return _1_8th * (1.0 + xi * xi_i[index]) *
                  (1.0 + eta * eta_i[index]) *
                  (1.0 + zeta * zeta_i[index]);
}

double HexahedronPWLFEValues::HexGradShape_x(unsigned int index,
                                             double xi,
                                             double eta,
                                             double zeta)
{
  return _1_8th * xi_i[index] *
                  (1.0 + eta * eta_i[index]) *
                  (1.0 + zeta * zeta_i[index]);
}

double HexahedronPWLFEValues::HexGradShape_y(unsigned int index,
                                             double xi,
                                             double eta,
                                             double zeta)
{
  return _1_8th * eta_i[index] *
         (1.0 + xi * xi_i[index]) *
         (1.0 + zeta * zeta_i[index]);
}

double HexahedronPWLFEValues::HexGradShape_z(unsigned int index,
                                             double xi,
                                             double eta,
                                             double zeta)
{
  return _1_8th * zeta_i[index] *
         (1.0 + xi * xi_i[index]) *
         (1.0 + eta * eta_i[index]);
}

