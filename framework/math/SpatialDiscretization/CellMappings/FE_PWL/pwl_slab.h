#ifndef PWL_SLAB_VALUES_H
#define PWL_SLAB_VALUES_H

#include "math/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"
#include "math/Quadratures/quadrature_line.h"
#include "mesh/Cell/cell.h"
#include "pwl_cellbase.h"
#include <array>


//###################################################################
namespace chi_math
{
  /**Object for handling slab shaped piecewise linear shape functions.*/
  class SlabMappingFE_PWL : public chi_math::CellMappingFE_PWL
  {
  private:
    chi_mesh::Vector3 v0_;
    uint64_t v0i_;
    uint64_t v1i_;
    std::array<chi_mesh::Normal,2> normals_;
    const QuadratureLine& volume_quadrature_;
    double h_;

  public:

    /**Constructor for a slab view.*/
    SlabMappingFE_PWL(const chi_mesh::Cell& slab_cell,
                      const chi_mesh::MeshContinuum& ref_grid,
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
    double ShapeValue(const int i, const chi_mesh::Vector3& xyz) const override;

    chi_mesh::Vector3 GradShapeValue(
      const int i,
      const chi_mesh::Vector3& xyz) const override;

    void ShapeValues(const chi_mesh::Vector3& xyz,
                     std::vector<double>& shape_values) const override;

    void GradShapeValues(
      const chi_mesh::Vector3& xyz,
      std::vector<chi_mesh::Vector3>& gradshape_values) const override;
  };
}
#endif
