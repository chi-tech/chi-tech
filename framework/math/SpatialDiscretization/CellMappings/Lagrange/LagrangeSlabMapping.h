#ifndef CHITECH_LAGRANGESLABMAPPING_H
#define CHITECH_LAGRANGESLABMAPPING_H

#include "math/SpatialDiscretization/CellMappings/LagrangeBaseMapping.h"

namespace chi_math::cell_mapping
{

/**Lagrange element mapping for a 1D Slab.
* \ingroup doc_CellMappings*/
class LagrangeSlabMapping : public LagrangeBaseMapping
{
public:
  LagrangeSlabMapping(const chi_mesh::MeshContinuum& grid,
                      const chi_mesh::Cell& cell,
                      const Quadrature& volume_quadrature,
                      const Quadrature& surface_quadrature,
                      CoordinateSystemType coordinate_system_type);

protected:
  double RefShape(uint32_t i, const Vec3& qpoint) const override;
  Vec3 RefGradShape(uint32_t i, const Vec3& qpoint) const override;

  MatDbl RefJacobian(const Vec3& qpoint) const override;

  std::pair<double, Vec3>
  RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                      const Vec3& qpoint_face) const override;

  Vec3 FaceToElementQPointConversion(size_t face_index,
                                     const Vec3& qpoint_face) const override;
};

} // namespace chi_math::cell_mapping

#endif // CHITECH_LAGRANGESLABMAPPING_H
