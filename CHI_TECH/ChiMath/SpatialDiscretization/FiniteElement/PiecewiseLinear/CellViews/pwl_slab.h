#ifndef PWL_SLAB_VALUES_H
#define PWL_SLAB_VALUES_H

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include <vector>
#include <ChiMesh/Cell/cell_slab.h>

//###################################################################
/**Object for handling slab shaped piecewise linear shape functions.*/
class SlabPWLFEView : public CellPWLFEValues
{
private:
  chi_mesh::Vector3 v0;
  int v0i;
  int v1i;
  std::array<chi_mesh::Normal,2> normals;
  const chi_math::QuadratureGaussLegendre& default_volume_quadrature;
  const chi_math::QuadratureGaussLegendre& arbitrary_volume_quadrature;
  double h;
public:

  /**Constructor for a slab view.*/
  SlabPWLFEView(const chi_mesh::CellSlab& slab_cell,
                const chi_mesh::MeshContinuumPtr& ref_grid,
                const chi_math::QuadratureGaussLegendre& minumum_volume_quadrature,
                const chi_math::QuadratureGaussLegendre& arb_volume_quadrature);

  void ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) override;
  void InitializeAllQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) override;

  void InitializeVolumeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data) override;

  void InitializeFaceQuadraturePointData(
    unsigned int face,
    chi_math::finite_element::FaceQuadraturePointData& faces_qp_data) override;

  //################################################## Define standard
  //                                                   slab linear shape
  //                                                   functions
  double SlabShape(int index, int qpoint_index, bool on_surface=false);
  double SlabGradShape(int index);


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
#endif
