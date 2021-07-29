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
    NO_FLAGS_SET = 0
  };

  inline SetupFlags
  operator|(const SetupFlags f1, const SetupFlags f2)
  {
    return static_cast<SetupFlags>(static_cast<unsigned int>(f1) |
                                   static_cast<unsigned int>(f2));
  }

  //#############################################
  /**Storage structure for unit integrals.*/
  class UnitIntegralData
  {
  public:
    typedef std::vector<double> VecDbl;
    typedef std::vector<VecDbl> MatDbl;
    typedef std::vector<chi_mesh::Vector3> VecVec3;
    typedef std::vector<VecVec3> MatVec3;

  private:
    MatDbl   m_IntV_gradShapeI_gradShapeJ;
    MatVec3  m_IntV_shapeI_gradshapeJ    ;
    MatDbl   m_IntV_shapeI_shapeJ        ;
    VecDbl   m_IntV_shapeI               ;
    VecVec3  m_IntV_gradshapeI           ;

    std::vector<MatDbl>  m_IntS_shapeI_shapeJ    ;
    std::vector<VecDbl>  m_IntS_shapeI           ;
    std::vector<MatVec3> m_IntS_shapeI_gradshapeJ;

    std::vector<std::vector<int>> m_face_dof_mappings;
    size_t m_num_nodes=0;

  public:
    void Initialize(MatDbl   in_IntV_gradShapeI_gradShapeJ,
                    MatVec3  in_IntV_shapeI_gradshapeJ,
                    MatDbl   in_IntV_shapeI_shapeJ,
                    VecDbl   in_IntV_shapeI,
                    VecVec3  in_IntV_gradshapeI,
                    std::vector<MatDbl>  in_IntS_shapeI_shapeJ,
                    std::vector<VecDbl>  in_IntS_shapeI,
                    std::vector<MatVec3> in_IntS_shapeI_gradshapeJ,
                    std::vector<std::vector<int>> in_face_dof_mappings,
                    size_t in_num_nodes);

    double IntV_gradShapeI_gradShapeJ(unsigned int i,
                                      unsigned int j) const;
    chi_mesh::Vector3 IntV_shapeI_gradshapeJ(unsigned int i,
                                             unsigned int j) const;
    double IntV_shapeI_shapeJ(unsigned int i,
                              unsigned int j) const;
    double IntV_shapeI(unsigned int i) const;
    chi_mesh::Vector3 IntV_gradshapeI(unsigned int i) const;

    double IntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j) const;

    double IntS_shapeI(unsigned int face, unsigned int i) const;

    chi_mesh::Vector3 IntS_shapeI_gradshapeJ(unsigned int face,
                                             unsigned int i,
                                             unsigned int j) const;
    int FaceDofMapping(size_t face, size_t face_node_index) const
    {
      auto& face_data = m_face_dof_mappings.at(face);
      return face_data.at(face_node_index);
    }
    size_t NumNodes() const
    {
      return m_num_nodes;
    }

    const MatDbl  & GetIntV_gradShapeI_gradShapeJ() const {return m_IntV_gradShapeI_gradShapeJ;}
    const MatVec3 & GetIntV_shapeI_gradshapeJ()     const {return m_IntV_shapeI_gradshapeJ    ;}
    const MatDbl  & GetIntV_shapeI_shapeJ()         const {return m_IntV_shapeI_shapeJ        ;}
    const VecDbl  & GetIntV_shapeI()                const {return m_IntV_shapeI               ;}
    const VecVec3 & GetIntV_gradshapeI()            const {return m_IntV_gradshapeI           ;}

    const std::vector<MatDbl>&  GetIntS_shapeI_shapeJ()    const  {return m_IntS_shapeI_shapeJ    ;}
    const std::vector<VecDbl>&  GetIntS_shapeI()           const  {return m_IntS_shapeI           ;}
    const std::vector<MatVec3>& GetIntS_shapeI_gradshapeJ()const  {return m_IntS_shapeI_gradshapeJ;}
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
    VecDbl                        m_JxW                     ; ///< qp index only
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
                        size_t num_nodes);
    const std::vector<unsigned int>&
      QuadraturePointIndices() const;
    chi_mesh::Vector3
      QPointXYZ(unsigned int qp) const;
    double
      ShapeValue(unsigned int i, unsigned int qp) const;
    chi_mesh::Vector3
      ShapeGrad(unsigned int i, unsigned int qp) const;
    double
      JxW(unsigned int qp) const;
    int
      FaceDofMapping(size_t face, size_t face_node_index) const;
    size_t
      NumNodes() const;
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
                        size_t num_nodes);
    chi_mesh::Vector3
      Normal(unsigned int qp) const;
  };
}

}

#endif //CHI_MATH_FINITE_ELEMENT_H
