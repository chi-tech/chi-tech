#ifndef PWL_CELL_VALUES_BASE_H
#define PWL_CELL_VALUES_BASE_H

#include <ChiMesh/chi_mesh.h>

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
/** Base class for all cell FE views.*/
class CellPWLFEValues
{
protected:
  chi_mesh::MeshContinuumPtr grid;
public:
  const int num_nodes;

  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;

  MatDbl   IntV_gradShapeI_gradShapeJ;
  MatVec3  IntV_shapeI_gradshapeJ;
  MatDbl   IntV_shapeI_shapeJ;
  VecDbl   IntV_shapeI;
  VecVec3  IntV_gradshapeI;

  std::vector<MatDbl>  IntS_shapeI_shapeJ;
  std::vector<VecDbl>  IntS_shapeI;
  std::vector<MatVec3> IntS_shapeI_gradshapeJ;

  std::vector<std::vector<int>> face_dof_mappings;

protected:
  bool precomputed = false;   ///< Are the integrals computed.

public:
  explicit CellPWLFEValues(int num_dofs,
                           chi_mesh::MeshContinuumPtr ref_grid) :
    grid(ref_grid),
    num_nodes(num_dofs)
  {}

  virtual void
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) {}

  virtual void
  InitializeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) {}

  virtual void PreComputeValues() {}

  /** Virtual function evaluation of the shape function. */
  virtual double ShapeValue(const int i, const chi_mesh::Vector3& xyz)
  {
    return 0.0;
  }

  /** Virtual function returning the all the shape function evaluations
   * at the point.*/
  virtual void ShapeValues(const chi_mesh::Vector3& xyz,
                           std::vector<double>& shape_values)
  {
    shape_values.resize(num_nodes, 0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual chi_mesh::Vector3 GradShapeValue(const int i,
                                           const chi_mesh::Vector3& xyz)
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual void GradShapeValues(const chi_mesh::Vector3& xyz,
                               std::vector<chi_mesh::Vector3>& gradshape_values)
  {
    gradshape_values.resize(num_nodes, chi_mesh::Vector3());
  }
  virtual ~CellPWLFEValues() = default;
};


#endif