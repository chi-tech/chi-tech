#ifndef SLAB_MAPPING_FE_PWL_CYLINDRICAL_H
#define SLAB_MAPPING_FE_PWL_CYLINDRICAL_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab.h"


/** Object for handling slab shaped 1D cells in axial-symmetric
 *  cylindrical coordinates. */
class SlabMappingFE_PWL_Cylindrical : public SlabMappingFE_PWL
{
//  Methods
public:
  SlabMappingFE_PWL_Cylindrical(
    const chi_mesh::Cell& slab_cell,
    const chi_mesh::MeshContinuumPtr& ref_grid,
    const chi_math::QuadratureLine& volume_quadrature)
  : SlabMappingFE_PWL(slab_cell, ref_grid, volume_quadrature)
  {}
private:
  double SpatialWeightFunction(const chi_mesh::Vector3& pt) const override
  { return pt[2]; }
public:
  void
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const override
  { ComputeWeightedUnitIntegrals(ui_data); }
};

#endif // SLAB_MAPPING_FE_PWL_CYLINDRICAL_H
