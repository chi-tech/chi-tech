#ifndef POLYGON_MAPPING_FE_PWL_CYLINDRICAL_H
#define POLYGON_MAPPING_FE_PWL_CYLINDRICAL_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon.h"


/** Object for handling polygon shaped 2D cells in axial-symmetric
 *  cylindrical coordinates. */
class PolygonMappingFE_PWL_Cylindrical : public PolygonMappingFE_PWL
{
//  Methods
public:
  PolygonMappingFE_PWL_Cylindrical(
    const chi_mesh::Cell& poly_cell,
    const chi_mesh::MeshContinuumPtr& ref_grid,
    const chi_math::QuadratureTriangle& volume_quadrature,
    const chi_math::QuadratureLine&     surface_quadrature)
  : PolygonMappingFE_PWL(poly_cell, ref_grid,
                         volume_quadrature, surface_quadrature)
  {}
private:
  double SpatialWeightFunction(const chi_mesh::Vector3& pt) const override
  { return pt[0]; }
public:
  void
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const override
  { ComputeWeightedUnitIntegrals(ui_data); }
};

#endif // POLYGON_MAPPING_FE_PWL_CYLINDRICAL_H
