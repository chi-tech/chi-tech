#ifndef CELL_FVDATA_BASE_H
#define CELL_FVDATA_BASE_H

#include <utility>

#include "ChiMath/SpatialDiscretization/CellMappings/cell_mapping_base.h"

#include "ChiMesh/chi_meshvector.h"

//######################################################### Class def
namespace chi_math
{

/**Base cell class for Finite Volume Method.*/
class CellFVValues : public CellMapping
{
protected:
  const chi_mesh::Vector3 m_cell_centroid;
public:
  double                volume=0.0;
  std::vector<double>   face_area; ///< Actually areas

  explicit CellFVValues(chi_mesh::MeshContinuumConstPtr grid,
                        const chi_mesh::Vector3& cc,
                        std::vector<std::vector<int>> face_node_mappings) :
    CellMapping(std::move(grid), 1, std::move(face_node_mappings)),
    m_cell_centroid(cc)
    {}

public:
  //02 Shapefuncs
  double ShapeValue(int i, const chi_mesh::Vector3& xyz)  override
  {
    return 1.0;
  }
  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values)  override
  {
    shape_values.assign(num_nodes, 0.0);
  }
  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz)  override
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }
  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values)
                        override
  {
    gradshape_values.assign(num_nodes, chi_mesh::Vector3(0,0,0));
  }
  std::vector<chi_mesh::Vector3> GetNodeLocations() const override
  {
    return {m_cell_centroid};
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