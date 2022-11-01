#ifndef CHITECH_CELL_MAPPING_BASE_H
#define CHITECH_CELL_MAPPING_BASE_H

#include <memory>
#include <utility>
#include <vector>

//################################################################### Fwd Decls.
namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<MeshContinuum const> MeshContinuumConstPtr;
  struct Vector3;
  class Cell;
}

namespace chi_math::finite_element
{
  class UnitIntegralData;
  class InternalQuadraturePointData;
  class FaceQuadraturePointData;
}

namespace chi_math
{
//################################################################### Class def
/**Base class for all cell mappings.*/
class CellMapping
{
protected:
  chi_mesh::MeshContinuumConstPtr grid;
  const size_t num_nodes;

  /** For each cell face, map from the face node index to the corresponding
     *  cell node index. More specifically, \p face_dof_mappings[f][fi], with
     *  \p fi the face node index of the face identified by face index \p f,
     *  contains the corresponding cell node index. */
  const std::vector<std::vector<int>> face_node_mappings;

  CellMapping(chi_mesh::MeshContinuumConstPtr  in_grid,
              size_t in_num_nodes,
              std::vector<std::vector<int>> in_face_node_mappings) :
              grid(std::move(in_grid)),
              num_nodes(in_num_nodes),
              face_node_mappings(std::move(in_face_node_mappings))
              {}

public:
  //00
  size_t NumNodes() const {return num_nodes;}

  //02 ShapeFuncs
  virtual
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const = 0;
  virtual
  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) const = 0;
  virtual
  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz) const = 0;
  virtual
  void GradShapeValues(
    const chi_mesh::Vector3& xyz,
    std::vector<chi_mesh::Vector3>& gradshape_values) const = 0;
  virtual
  std::vector<chi_mesh::Vector3> GetNodeLocations() const = 0;

  //03 Quadrature
  /** Compute unit integrals. */
  virtual void
  ComputeUnitIntegrals(finite_element::UnitIntegralData& ui_data) const;

  /** Initialize volume quadrature point data and
   *  surface quadrature point data for all faces. */
  void
  InitializeAllQuadraturePointData(
    finite_element::InternalQuadraturePointData& internal_data,
    std::vector<finite_element::FaceQuadraturePointData>& faces_qp_data) const;

  /** Initialize volume quadrature point data. */
  virtual void
  InitializeVolumeQuadraturePointData(
    finite_element::InternalQuadraturePointData& internal_data) const = 0;

  /** Initialize surface quadrature point data for face index \p face. */
  virtual void
  InitializeFaceQuadraturePointData(
    unsigned int face,
    finite_element::FaceQuadraturePointData& faces_qp_data) const = 0;

public:
  virtual ~CellMapping() = default;
};
}//namespace chi_math

#endif //CHITECH_CELL_MAPPING_BASE_H
