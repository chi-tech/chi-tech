#ifndef _product_quadrature_h
#define _product_quadrature_h

#include <vector>
#include "../../ChiMesh/chi_mesh.h"

#define GAUSS_LEGENDRE           1
#define GAUSS_CHEBYSHEV          2
#define GAUSS_LEGENDRE_LEGENDRE  3
#define GAUSS_LEGENDRE_CHEBYSHEV 4
#define CUSTOM_QUADRATURE        5

#include "angular_quadrature_base.h"

namespace chi_math
{
  class ProductQuadrature;

}

//######################################################### Class def
/** Class for product quadratures*/
class chi_math::ProductQuadrature : public chi_math::AngularQuadrature
{
public:
  std::vector<double>           polar_ang;
  std::vector<double>           azimu_ang;

public:
  //product_quadrature.cc
       ProductQuadrature() :
         AngularQuadrature(chi_math::AngularQuadratureType::ProductQuadrature)
       {}
  void InitializeWithGL(int Np, bool verbose=false);
  void InitializeWithGLL(int Na, int Np, bool verbose=false);
  void InitializeWithGLC(int Na, int Np, bool verbose=false);
  void InitializeWithCustom(std::vector<double>& azimuthal,
                            std::vector<double>& polar,
                            std::vector<double>& in_weights,
                            bool verbose=false) override;
  int  GetAngleNum(int polar_angle_index, int azimu_ang_index);

};




#endif