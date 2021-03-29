#ifndef CELL_MAPPING_FE_PWL_BASE_H
#define CELL_MAPPING_FE_PWL_BASE_H

#include "ChiMesh/chi_mesh.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
/** Base class for all cell FE views.*/
class CellMappingFE_PWL
{
protected:
  chi_mesh::MeshContinuumPtr grid;
public:
  const int num_nodes;

  typedef std::vector<double> VecDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;

  std::vector<std::vector<int>> face_dof_mappings;

public:
  explicit CellMappingFE_PWL(int num_dofs,
                             chi_mesh::MeshContinuumPtr ref_grid) :
    grid(std::move(ref_grid)),
    num_nodes(num_dofs)
  {}

  virtual void
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const;

  virtual void
  InitializeAllQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) const;

  virtual void
  InitializeVolumeQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data) const = 0;

  virtual void
  InitializeFaceQuadraturePointData(
    unsigned int face,
    chi_math::finite_element::FaceQuadraturePointData& faces_qp_data) const = 0;

public:
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
  virtual ~CellMappingFE_PWL() = default;

protected:
  /** Spatial weight function. See also ComputeWeightedUnitIntegrals. */
  virtual double SpatialWeightFunction(const chi_mesh::Vector3& pt) const
  { return static_cast<double>(0); }
  /** Compute spatially-weighted unit integrals from quadrature point data.
   *  Spatial weighting is given by SpatialWeightFunction evaluated at
   *  the real quadrature points. */
  void
  ComputeWeightedUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const;
};


#endif