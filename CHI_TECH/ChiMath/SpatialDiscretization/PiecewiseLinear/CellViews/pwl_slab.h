#ifndef PWL_SLAB_VALUES_H
#define PWL_SLAB_VALUES_H

#include "../pwl.h"
#include <vector>
#include <ChiMesh/Cell/cell_slab.h>

//###################################################################
/**Object for handling slab shaped piecewise linear shape functions.*/
class SlabPWLFEView : public CellPWLFEValues
{
private:
  chi_mesh::MeshContinuumPtr grid;
  int v0i;
  int v1i;
  chi_math::QuadratureGaussLegendre default_volume_quadrature;
public:
  double h;
public:

  /**Constructor for a slab view.*/
  SlabPWLFEView(chi_mesh::CellSlab *slab_cell,
                chi_mesh::MeshContinuumPtr& in_grid) :
    CellPWLFEValues(2),
    default_volume_quadrature(chi_math::QuadratureGaussLegendre(chi_math::QuadratureOrder::SECOND))
  {
    grid = in_grid;
    v0i = slab_cell->vertex_ids[0];
    v1i = slab_cell->vertex_ids[1];
    chi_mesh::Vertex v0 = *grid->vertices[v0i];
    chi_mesh::Vertex v1 = *grid->vertices[v1i];

    chi_mesh::Vector3 v01 = v1 - v0;
    h = v01.Norm();

    face_dof_mappings.emplace_back(1,0);
    face_dof_mappings.emplace_back(1,1);

  }

  void ComputeUnitIntegrals();

  //################################################## Define standard
  //                                                   slab linear shape
  //                                                   functions
  double SlabShape(int index, int qpoint_index, bool on_surface=false);

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

  void PreComputeValues() override;
};
#endif
