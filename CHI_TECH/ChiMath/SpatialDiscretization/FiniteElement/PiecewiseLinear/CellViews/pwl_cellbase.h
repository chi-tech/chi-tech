#ifndef PWL_CELL_VALUES_BASE_H
#define PWL_CELL_VALUES_BASE_H

#include <ChiMesh/chi_mesh.h>

#define PWL_CELL_THROW_QP_UNINIT throw std::invalid_argument(\
                                       "InternalQuadraturePointData called "\
                                       "without being initialized.")

//###################################################################
/** Base class for all cell FE views.*/
class CellPWLFEValues
{
protected:
  chi_mesh::MeshContinuumPtr grid;
public:
  const int dofs;

  std::vector<std::vector<double>>              IntV_gradShapeI_gradShapeJ;
  std::vector<std::vector<chi_mesh::Vector3>>   IntV_shapeI_gradshapeJ;
  std::vector<std::vector<double>>              IntV_shapeI_shapeJ;
  std::vector<double>                           IntV_shapeI;

  std::vector<std::vector<std::vector<double>>> IntS_shapeI_shapeJ;
  std::vector<std::vector<double>>              IntS_shapeI;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> IntS_shapeI_gradshapeJ;

  std::vector<std::vector<int>> face_dof_mappings;

protected:
  bool precomputed = false;   ///< Are the integrals computed.

public:
  typedef std::vector<double> VecDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;

  class InternalQuadraturePointData
  {
    friend class SlabPWLFEView;
    friend class PolygonPWLFEValues;
    friend class PolyhedronPWLFEValues;
  public:
    std::vector<unsigned int> quadrature_point_indices; ///< qp only
  protected:
    std::vector<VecDbl>       m_shape_value;              ///< Node i, then qp
    std::vector<VecVec3>      m_shape_grad;               ///< Node i, then qp
    VecDbl                    m_JxW;                      ///< Node i, then qp
    bool                      initialized=false;

  public:
    double shape_value(unsigned int i, unsigned int qp) const
    {
      if (not initialized) PWL_CELL_THROW_QP_UNINIT;
      auto& qp_data = m_shape_value.at(i);
      return qp_data.at(qp);
    }
    chi_mesh::Vector3 shape_grad(unsigned int i, unsigned int qp) const
    {
      if (not initialized) PWL_CELL_THROW_QP_UNINIT;
      auto& qp_data = m_shape_grad.at(i);
      return qp_data.at(qp);
    }
    double JxW(unsigned int qp) const
    {
      if (not initialized) PWL_CELL_THROW_QP_UNINIT;
      return m_JxW.at(qp);
    }
  }internal_qp_data;

  class FaceQuadraturePointData : InternalQuadraturePointData
  {
    friend class SlabPWLFEView;
    friend class PolygonPWLFEValues;
    friend class PolyhedronPWLFEValues;
  public:
    VecVec3                   m_normals;                  ///< node i, then qp
    double normal(unsigned int qp) const
    {
      if (not initialized) PWL_CELL_THROW_QP_UNINIT;
      return m_JxW.at(qp);
    }
  };
  std::vector<FaceQuadraturePointData> surface_qp_data;

public:
  explicit CellPWLFEValues(int num_dofs,
                           chi_mesh::MeshContinuumPtr ref_grid) :
                           grid(ref_grid),
                           dofs(num_dofs)
  {}

  virtual ~CellPWLFEValues() = default;

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
    shape_values.resize(dofs,0.0);
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
    gradshape_values.resize(dofs,chi_mesh::Vector3());
  }

  virtual void PreComputeValues() {}
};


#endif