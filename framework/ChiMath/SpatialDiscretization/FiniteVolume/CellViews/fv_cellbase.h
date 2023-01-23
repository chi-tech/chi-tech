#ifndef CELL_FVDATA_BASE_H
#define CELL_FVDATA_BASE_H

#include <utility>

#include "ChiMath/SpatialDiscretization/CellMappings/cell_mapping_base.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//######################################################### Class def
namespace chi_math
{

/**Base cell class for Finite Volume Method.*/
class CellFVValues : public CellMapping
{
public:
  explicit CellFVValues(chi_mesh::MeshContinuumConstPtr grid,
                        const chi_mesh::Cell& cell,
                        const chi_mesh::Vector3& cc,
                        std::vector<std::vector<int>> face_node_mappings) :
    CellMapping(std::move(grid), cell, 1, std::move(face_node_mappings),
                &CellMapping::ComputeCellVolumeAndAreas)
    {}

public:
  //02 Shapefuncs
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override
  {
    if (m_grid_ptr->CheckPointInsideCell(m_cell, xyz))
      return 1.0;
    else
      return 0.0;
  }
  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) const override
  {
    if (m_grid_ptr->CheckPointInsideCell(m_cell, xyz))
      shape_values.assign(m_num_nodes, 1.0);
    else
      shape_values.assign(m_num_nodes, 0.0);
  }
  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz) const override
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }
  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values)
                       const override
  {
    gradshape_values.assign(m_num_nodes, chi_mesh::Vector3(0, 0, 0));
  }
  std::vector<chi_mesh::Vector3> GetNodeLocations() const override
  {
    return {m_cell.centroid};
  }

  //03 Quadrature
  void InitializeVolumeQuadraturePointData(
    finite_element::InternalQuadraturePointData& internal_data) const override;


  void InitializeFaceQuadraturePointData(
    unsigned int face,
    finite_element::FaceQuadraturePointData& faces_qp_data) const override;
};

}//namespace chi_math

#endif