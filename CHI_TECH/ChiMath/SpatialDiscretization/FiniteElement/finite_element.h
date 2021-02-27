#ifndef CHI_MATH_FINITE_ELEMENT_H
#define CHI_MATH_FINITE_ELEMENT_H

#include "ChiMath/chi_math.h"

namespace chi_math
{

namespace finite_element
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;

  //#############################################
  enum SetupFlags : int
  {
    NO_FLAGS_SET           = 0,
    COMPUTE_UNIT_INTEGRALS = (1 << 0),
    INIT_QP_DATA           = (1 << 1)
  };

  //#############################################
  /**Storage structure for unit integrals.*/
  class UnitIntegralData
  {
  public:
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
    size_t num_nodes=0;
  };

  //#############################################
  /**Stored relevant quadrature point information
   * for volumetric integrals.*/
  class InternalQuadraturePointData
  {
  protected:
    static void THROW_QP_UNINIT()
    {throw std::invalid_argument("InternalQuadraturePointData called "
                                 "without being initialized. Set flag"
                                 " INIT_QP_DATA.");}
  protected:
    std::vector<unsigned int>     m_quadrature_point_indices; ///< qp index only
    VecVec3                       m_qpoints_xyz             ; ///< qp index only
    std::vector<VecDbl>           m_shape_value             ; ///< Node i, then qp
    std::vector<VecVec3>          m_shape_grad              ; ///< Node i, then qp
    VecDbl                        m_JxW                     ; ///< Node i, then qp
    std::vector<std::vector<int>> m_face_dof_mappings       ; ///< Face f,then fi
    size_t                        m_num_nodes                =0;

    bool                          m_initialized=false;

  public:
    void InitializeData(std::vector<unsigned int> quadrature_point_indices,
                        VecVec3                   qpoints_xyz,
                        std::vector<VecDbl>       shape_value,
                        std::vector<VecVec3>      shape_grad,
                        VecDbl                    JxW,
                        std::vector<std::vector<int>> face_dof_mappings,
                        size_t num_nodes)
    {
      m_quadrature_point_indices= std::move(quadrature_point_indices);
      m_qpoints_xyz             = std::move(qpoints_xyz             );
      m_shape_value             = std::move(shape_value             );
      m_shape_grad              = std::move(shape_grad              );
      m_JxW                     = std::move(JxW                     );
      m_face_dof_mappings       = std::move(face_dof_mappings       );
      m_num_nodes               = num_nodes               ;
      m_initialized = true;
    }
    const std::vector<unsigned int>& quadrature_point_indices() const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      return m_quadrature_point_indices;
    }
    chi_mesh::Vector3 qpoint_xyz(unsigned int qp) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      return m_qpoints_xyz.at(qp);
    }
    double shape_value(unsigned int i, unsigned int qp) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      auto& qp_data = m_shape_value.at(i);
      return qp_data.at(qp);
    }
    chi_mesh::Vector3 shape_grad(unsigned int i, unsigned int qp) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      auto& qp_data = m_shape_grad.at(i);
      return qp_data.at(qp);
    }
    double JxW(unsigned int qp) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      return m_JxW.at(qp);
    }
    int face_dof_mapping(size_t face, size_t face_node_index) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      auto& face_data = m_face_dof_mappings.at(face);
      return face_data.at(face_node_index);
    }
    size_t num_nodes() const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      return m_num_nodes;
    }
  };

  //#############################################
  /**Stores relevant quadrature point information
   * for surface integrals.*/
  class FaceQuadraturePointData : public InternalQuadraturePointData
  {
  protected:
    VecVec3                   m_normals;                  ///< node i, then qp
  public:
    void InitializeData(std::vector<unsigned int>     quadrature_point_indices,
                        VecVec3                       qpoints_xyz,
                        std::vector<VecDbl>           shape_value,
                        std::vector<VecVec3>          shape_grad,
                        VecDbl                        JxW,
                        VecVec3                       normals,
                        std::vector<std::vector<int>> face_dof_mappings,
                        size_t num_nodes)
    {
      m_quadrature_point_indices= std::move(quadrature_point_indices);
      m_qpoints_xyz             = std::move(qpoints_xyz             );
      m_shape_value             = std::move(shape_value             );
      m_shape_grad              = std::move(shape_grad              );
      m_JxW                     = std::move(JxW                     );
      m_normals                 = std::move(normals                 );
      m_face_dof_mappings       = std::move(face_dof_mappings       );
      m_num_nodes               = num_nodes               ;
      m_initialized = true;
    }
    double normal(unsigned int qp) const
    {
      if (not m_initialized) THROW_QP_UNINIT();
      return m_JxW.at(qp);
    }
  };
}

}

#endif //CHI_MATH_FINITE_ELEMENT_H