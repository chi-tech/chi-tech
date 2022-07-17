#ifndef PWL_SLAB_VALUES_H
#define PWL_SLAB_VALUES_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"
#include "ChiMath/Quadratures/quadrature_line.h"
#include "ChiMesh/Cell/cell.h"
#include "pwl_cellbase.h"
#include <array>


//###################################################################
namespace chi_math
{
  /**Object for handling slab shaped piecewise linear shape functions.*/
  class SlabMappingFE_PWL : public chi_math::CellMappingFE_PWL
  {
  private:
    chi_mesh::Vector3 v0;
    int v0i;
    int v1i;
    std::array<chi_mesh::Normal,2> normals;
    const QuadratureLine& volume_quadrature;
    double h;
  public:

    /**Constructor for a slab view.*/
    SlabMappingFE_PWL(const chi_mesh::Cell& slab_cell,
                      const chi_mesh::MeshContinuumPtr& ref_grid,
                      const QuadratureLine& volume_quadrature);

    void ComputeUnitIntegrals(
      finite_element::UnitIntegralData& ui_data) const override;

    void InitializeVolumeQuadraturePointData(
      finite_element::InternalQuadraturePointData& internal_data) const override;

    void InitializeFaceQuadraturePointData(
      unsigned int face,
      finite_element::FaceQuadraturePointData& faces_qp_data) const override;

    //################################################## Define standard
    //                                                   slab linear shape
    //                                                   functions
    double SlabShape(int index,
                     const chi_mesh::Vector3& qpoint,
                     bool on_surface=false,
                     const int edge=0) const;
    double SlabGradShape(int index) const;


    //############################################### Actual shape functions
    //                                                as function of cartesian
    //                                                coordinates
  public:
    double ShapeValue(const int i, const chi_mesh::Vector3& xyz) override;
    chi_mesh::Vector3 GradShapeValue(const int i, const chi_mesh::Vector3& xyz) override;

    void ShapeValues(const chi_mesh::Vector3& xyz,
                     std::vector<double>& shape_values) override;
    void GradShapeValues(const chi_mesh::Vector3& xyz,
                         std::vector<chi_mesh::Vector3>& gradshape_values) override;


  };
}
#endif
