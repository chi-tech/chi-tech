#ifndef _angular_quadrature_base_h
#define _angular_quadrature_base_h

#include <vector>

#include "ChiMesh/chi_mesh.h"

namespace chi_math
{
  struct QuadraturePointPhiTheta;

  enum class AngularQuadratureType
  {
    Arbitrary         = 1,
    ProductQuadrature = 2,
    SLDFESQ           = 3
  };
  class AngularQuadrature;
  class AngularQuadratureCustom;
}

/**Simple structure to add names to the angle components.*/
struct chi_math::QuadraturePointPhiTheta
{
  double phi=0.0;
  double theta=0.0;
  QuadraturePointPhiTheta(const double phi, const double theta)
  : phi(phi), theta(theta) {}
};

//################################################################### Class def
/**Base class for angular quadratures.*/
class chi_math::AngularQuadrature
{
public:
  const chi_math::AngularQuadratureType type;
public:
  std::vector<chi_math::QuadraturePointPhiTheta> abscissae;
  std::vector<double>                            weights;
  std::vector<chi_mesh::Vector3>                 omegas;

  struct HarmonicIndices
  {
    unsigned int ell=0;
    int          m=0;

    HarmonicIndices()=default;
    HarmonicIndices(unsigned int in_ell, int in_m) : ell(in_ell),m(in_m){}

    bool operator==(const HarmonicIndices& other) const
    {
      return (ell == other.ell and m == other.m);
    }
  };

protected:
  std::vector<std::vector<double>> d2m_op;
  std::vector<std::vector<double>> m2d_op;
  std::vector<HarmonicIndices>     m_to_ell_em_map;
  bool                             d2m_op_built = false;
  bool                             m2d_op_built = false;

public:
  AngularQuadrature() :
  type(chi_math::AngularQuadratureType::Arbitrary)
  {}

  explicit
  AngularQuadrature(chi_math::AngularQuadratureType in_type) :
    type(in_type)
  {}

  virtual ~AngularQuadrature() = default;

  virtual void OptimizeForPolarSymmetry(double normalization);

protected:
  virtual void MakeHarmonicIndices(unsigned int scattering_order, int dimension);

public:
  virtual void BuildDiscreteToMomentOperator(unsigned int scattering_order,
                                             int dimension);
  virtual void BuildMomentToDiscreteOperator(unsigned int scattering_order,
                                             int dimension);

  std::vector<std::vector<double>> const&
  GetDiscreteToMomentOperator() const;

  std::vector<std::vector<double>> const&
  GetMomentToDiscreteOperator() const;

  const std::vector<HarmonicIndices>&
  GetMomentToHarmonicsIndexMap() const;
};

class chi_math::AngularQuadratureCustom : public chi_math::AngularQuadrature
{
public:
  AngularQuadratureCustom(std::vector<double>& azimuthal,
                          std::vector<double>& polar,
                          std::vector<double>& in_weights, bool verbose);
};

#endif
