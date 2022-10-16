#ifndef CELL_FVDATA_BASE_H
#define CELL_FVDATA_BASE_H

#include "ChiMath/SpatialDiscretization/CellMappings/cell_mapping_base.h"

//######################################################### Class def
namespace chi_math
{

/**Base cell class for Finite Volume Method.*/
class CellFVValues : public CellMapping
{
public:
  double                volume=0.0;
  std::vector<double>   face_area; ///< Actually areas

  explicit CellFVValues(chi_mesh::MeshContinuumPtr grid, size_t num_nodes) :
    CellMapping(std::move(grid), num_nodes)
    {}

public:
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
};

}//namespace chi_math

#endif