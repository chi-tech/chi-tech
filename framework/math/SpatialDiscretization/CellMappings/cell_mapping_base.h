#ifndef CHITECH_CELL_MAPPING_BASE_H
#define CHITECH_CELL_MAPPING_BASE_H

#include <memory>
#include <utility>
#include <vector>
#include <functional>

//################################################################### Fwd Decls.
namespace chi_mesh
{
  class MeshContinuum;
  typedef std::shared_ptr<const MeshContinuum> MeshContinuumConstPtr;
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
  /**This function gets called to compute the cell-volume and
   * face-areas. If simple linear cells are used then the
   * default CellMapping::ComputeCellVolumeAndAreas can be
   * used as this function. Otherwise (i.e. for higher order
   * elements, the child-class should
   * bind a different function to this.*/
  typedef std::function<void(const chi_mesh::MeshContinuum&,
                             const chi_mesh::Cell&,
                             double&,
                             std::vector<double>&)> VandAFunction;
protected:
  const chi_mesh::MeshContinuum& ref_grid_;
  const chi_mesh::Cell& cell_;
  const size_t num_nodes_;
  double volume_ = 0.0;
  std::vector<double> areas_;

  /** For each cell face, map from the face node index to the corresponding
     *  cell node index. More specifically, \p face_dof_mappings[f][fi], with
     *  \p fi the face node index of the face identified by face index \p f,
     *  contains the corresponding cell node index. */
  const std::vector<std::vector<int>> face_node_mappings_;

  CellMapping(const chi_mesh::MeshContinuum& in_grid,
              const chi_mesh::Cell& in_cell,
              size_t in_num_nodes,
              std::vector<std::vector<int>> in_face_node_mappings,
              const VandAFunction& volume_area_function);

public:
  //00
  size_t NumNodes() const { return num_nodes_; }
  size_t NumFaceNodes(size_t face_index) const
  {
    return face_node_mappings_.at(face_index).size();
  }

  static void ComputeCellVolumeAndAreas(
    const chi_mesh::MeshContinuum& grid,
    const chi_mesh::Cell& cell,
    double& volume,
    std::vector<double>& areas);

  double CellVolume() const {return volume_;}
  double FaceArea(size_t face_index) const {return areas_[face_index];}

  int MapFaceNode(size_t face_index, size_t face_node_index) const;

  //02 ShapeFuncs
  virtual double ShapeValue(int i, const chi_mesh::Vector3& xyz) const = 0;

  virtual void ShapeValues(const chi_mesh::Vector3& xyz,
                           std::vector<double>& shape_values) const = 0;
  virtual chi_mesh::Vector3
  GradShapeValue(int i, const chi_mesh::Vector3& xyz) const = 0;

  virtual void GradShapeValues(
      const chi_mesh::Vector3& xyz,
      std::vector<chi_mesh::Vector3>& gradshape_values) const = 0;

  virtual std::vector<chi_mesh::Vector3> GetNodeLocations() const = 0;

  //03 Quadrature
  /** Compute unit integrals. */
  virtual void ComputeUnitIntegrals(
      finite_element::UnitIntegralData& ui_data) const;

  /** Initialize volume quadrature point data and
   *  surface quadrature point data for all faces. */
  void InitializeAllQuadraturePointData(
    finite_element::InternalQuadraturePointData& internal_data,
    std::vector<finite_element::FaceQuadraturePointData>& faces_qp_data) const;

  /** Initialize volume quadrature point data. */
  virtual void
  InitializeVolumeQuadraturePointData(
    finite_element::InternalQuadraturePointData& internal_data) const = 0;

  /** Initialize surface quadrature point data for face index \p face. */
  virtual void InitializeFaceQuadraturePointData(
    unsigned int face,
    finite_element::FaceQuadraturePointData& faces_qp_data) const = 0;

  finite_element::InternalQuadraturePointData
  MakeVolumeQuadraturePointData() const;

  finite_element::FaceQuadraturePointData
  MakeFaceQuadraturePointData(size_t face_index) const;

public:
  virtual ~CellMapping() = default;
};
}//namespace chi_math

#endif //CHITECH_CELL_MAPPING_BASE_H
