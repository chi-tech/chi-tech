#ifndef _product_quadrature_h
#define _product_quadrature_h

#include <vector>

#include "angular_quadrature_base.h"

namespace chi_math
{
  enum class ProductQuadratureType
  {
    UNKNOWN                  = 0,
    GAUSS_LEGENDRE           = 1,
    GAUSS_CHEBYSHEV          = 2,
    GAUSS_LEGENDRE_LEGENDRE  = 3,
    GAUSS_LEGENDRE_CHEBYSHEV = 4,
    CUSTOM_QUADRATURE        = 5,
  };

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
  ProductQuadrature() :
    AngularQuadrature(chi_math::AngularQuadratureType::ProductQuadrature)
  {}

  virtual ~ProductQuadrature()
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
