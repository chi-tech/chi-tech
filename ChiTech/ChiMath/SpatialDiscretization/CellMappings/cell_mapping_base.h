#ifndef CHITECH_CELL_MAPPING_BASE_H
#define CHITECH_CELL_MAPPING_BASE_H

#include <memory>
#include <utility>

namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
}

namespace chi_math
{

class CellMapping
{
protected:
  chi_mesh::MeshContinuumPtr grid;
  const size_t num_nodes;

  CellMapping(chi_mesh::MeshContinuumPtr  in_grid,
              size_t in_num_nodes) :
              grid(std::move(in_grid)),
              num_nodes(in_num_nodes) {}

public:
  //02 ShapeFuncs
  virtual
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) = 0;
  virtual
  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) = 0;
  virtual
  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz) = 0;
  virtual
  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) = 0;

public:
  virtual ~CellMapping() = default;
};
}//namespace chi_math

#endif //CHITECH_CELL_MAPPING_BASE_H
