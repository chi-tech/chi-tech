#ifndef CELL_MAPPING_FE_PWL_BASE_H
#define CELL_MAPPING_FE_PWL_BASE_H

#include "ChiMesh/chi_mesh.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
namespace chi_math
{
  /** Base class for all cell FE views.*/
  class CellMappingFE_PWL
  {
  protected:
    /** Mesh. */
    chi_mesh::MeshContinuumPtr grid;
  public:
    typedef std::vector<double> VecDbl;
    typedef std::vector<chi_mesh::Vector3> VecVec3;
    /** Number of nodes of the cell. */
    const int num_nodes;
    /** For each cell face, map from the face node index to the corresponding
     *  cell node index. More specifically, \p face_dof_mappings[f][fi], with
     *  \p fi the face node index of the face identified by face index \p f,
     *  contains the corresponding cell node index. */
    std::vector<std::vector<int>> face_dof_mappings;

  public:
    /** Constructor. */
    explicit CellMappingFE_PWL(int num_dofs,
                               chi_mesh::MeshContinuumPtr ref_grid) :
      grid(std::move(ref_grid)),
      num_nodes(num_dofs)
    {}

    /** Compute unit integrals. */
    virtual void
    ComputeUnitIntegrals(finite_element::UnitIntegralData& ui_data) const;

    /** Initialize volume quadrature point data and
     *  surface quadrature point data for all faces. */
    virtual void
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

    /** Destructor. */
    virtual ~CellMappingFE_PWL() = default;

  protected:
    /** Spatial weight function. See also ComputeWeightedUnitIntegrals. */
    virtual double SpatialWeightFunction(const chi_mesh::Vector3& pt) const
    { return 0.0; }
    /** Compute spatially-weighted unit integrals from quadrature point data.
     *  Spatial weighting is given by SpatialWeightFunction evaluated at
     *  the real quadrature points. */
    void
    ComputeWeightedUnitIntegrals(finite_element::UnitIntegralData& ui_data) const;
  };
}


#endif