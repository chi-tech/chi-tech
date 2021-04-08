#ifndef _product_quadrature_h
#define _product_quadrature_h

#include <map>
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
protected:
  /** Linear indices of ordered directions mapped to polar level. */
  std::map<unsigned int, std::vector<unsigned int>> map_directions;

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
  /**Obtains the abscissae index given the indices of the
   * polar angle index and the azimuthal angle index.*/
  unsigned int GetAngleNum(const unsigned int polar_angle_index,
                           const unsigned int azimu_angle_index) const
  { return map_directions.at(polar_angle_index)[azimu_angle_index]; }
  /** Return constant reference to map_directions. */
  const std::map<unsigned int,
                 std::vector<unsigned int>>& GetDirectionMap() const
  { return map_directions; }
};

#endif
