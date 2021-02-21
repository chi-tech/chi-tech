#ifndef CHI_MATH_FINITE_ELEMENT_H
#define CHI_MATH_FINITE_ELEMENT_H

#include "ChiMath/chi_math.h"

#define PWL_CELL_THROW_QP_UNINIT throw std::invalid_argument(\
                                       "InternalQuadraturePointData called "\
                                       "without being initialized.")

namespace chi_math
{

namespace finite_element
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;

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
  public:
    std::vector<unsigned int> quadrature_point_indices; ///< qp only
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
  };

  //#############################################
  /**Stores relevant quadrature point information
   * for surface integrals.*/
  class FaceQuadraturePointData : public InternalQuadraturePointData
  {
  public:
    VecVec3                   m_normals;                  ///< node i, then qp
    double normal(unsigned int qp) const
    {
      if (not initialized) PWL_CELL_THROW_QP_UNINIT;
      return m_JxW.at(qp);
    }
  };
}

}

#endif //CHI_MATH_FINITE_ELEMENT_H