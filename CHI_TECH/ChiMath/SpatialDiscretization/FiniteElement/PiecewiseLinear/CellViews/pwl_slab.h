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
  chi_math::QuadratureGaussLegendre default_volume_quadrature;
  chi_math::QuadratureGaussLegendre arbitrary_volume_quadrature;
  double h;
public:

  /**Constructor for a slab view.*/
  SlabPWLFEView(chi_mesh::CellSlab *slab_cell,
                chi_mesh::MeshContinuumPtr ref_grid,
                chi_math::QuadratureGaussLegendre& minumum_volume_quadrature,
                chi_math::QuadratureGaussLegendre& arb_volume_quadrature) :
    CellPWLFEValues(2,ref_grid),
    default_volume_quadrature(minumum_volume_quadrature),
    arbitrary_volume_quadrature(arb_volume_quadrature)
  {
    grid = ref_grid;
    v0i = slab_cell->vertex_ids[0];
    v1i = slab_cell->vertex_ids[1];
                     v0 = *grid->vertices[v0i];
    chi_mesh::Vertex v1 = *grid->vertices[v1i];

    chi_mesh::Vector3 v01 = v1 - v0;
    h = v01.Norm();

    face_dof_mappings.emplace_back(1,0);
    face_dof_mappings.emplace_back(1,1);

    normals[0] = slab_cell->faces[0].normal;
    normals[1] = slab_cell->faces[1].normal;
  }

  void ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) override;
  void InitializeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) override;

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
