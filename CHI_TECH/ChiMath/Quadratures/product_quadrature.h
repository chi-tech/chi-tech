#ifndef _product_quadrature_h
#define _product_quadrature_h

#include <vector>
#include "../../ChiMesh/chi_mesh.h"

#define GAUSS_LEGENDRE           1
#define GAUSS_CHEBYSHEV          2
#define GAUSS_LEGENDRE_LEGENDRE  3
#define GAUSS_LEGENDRE_CHEBYSHEV 4

namespace chi_math
{
  struct QuadraturePointPhiTheta;
  class ProductQuadrature;
}

struct chi_math::QuadraturePointPhiTheta
{
  double phi;
  double theta;
};

//######################################################### Class def
/** Parent class for product quadratures*/
class chi_math::ProductQuadrature
{
public:
  std::vector<chi_math::QuadraturePointPhiTheta*> abscissae;
  std::vector<double>           polar_ang;
  std::vector<double>           azimu_ang;
  std::vector<double>           weights;
  std::vector<chi_mesh::Vector*>          omegas;

public:
  //product_quadrature.cc
  void InitializeWithGL(int Np, bool verbose=false);
  void InitializeWithGLL(int Na, int Np, bool verbose=false);
  void InitializeWithGLC(int Na, int Np, bool verbose=false);
  int  GetAngleNum(int polar_angle_index, int azimu_ang_index)
  {
    return azimu_ang_index*polar_ang.size() +
           polar_angle_index;
  }

};




#endif