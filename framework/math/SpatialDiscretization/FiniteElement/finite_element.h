#ifndef CHI_MATH_FINITE_ELEMENT_H
#define CHI_MATH_FINITE_ELEMENT_H

#include "math/chi_math.h"

namespace chi_math::finite_element
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;

  //#############################################
  enum SetupFlags : int
  {
    NO_FLAGS_SET           = 0,
    COMPUTE_CELL_MAPPINGS  = (1 << 0),
    COMPUTE_UNIT_INTEGRALS = (1 << 1),
    COMPUTE_QP_DATA        = (1 << 2)
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
    MatDbl   IntV_gradShapeI_gradShapeJ_;
    MatVec3  IntV_shapeI_gradshapeJ_    ;
    MatDbl   IntV_shapeI_shapeJ_        ;
    VecDbl   IntV_shapeI_               ;
    VecVec3  IntV_gradshapeI_           ;

    std::vector<MatDbl>  IntS_shapeI_shapeJ_    ;
    std::vector<VecDbl>  IntS_shapeI_           ;
    std::vector<MatVec3> IntS_shapeI_gradshapeJ_;

    std::vector<std::vector<int>> face_dof_mappings_;
    size_t num_nodes_ = 0;

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

    void Reset();

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
      auto& face_data = face_dof_mappings_.at(face);
      return face_data.at(face_node_index);
    }

    size_t NumNodes() const { return num_nodes_; }

    const MatDbl  & GetIntV_gradShapeI_gradShapeJ() const {return IntV_gradShapeI_gradShapeJ_;}
    const MatVec3 & GetIntV_shapeI_gradshapeJ()     const {return IntV_shapeI_gradshapeJ_    ;}
    const MatDbl  & GetIntV_shapeI_shapeJ()         const {return IntV_shapeI_shapeJ_        ;}
    const VecDbl  & GetIntV_shapeI()                const {return IntV_shapeI_               ;}
    const VecVec3 & GetIntV_gradshapeI()            const {return IntV_gradshapeI_           ;}

    const std::vector<MatDbl>&  GetIntS_shapeI_shapeJ()    const  {return IntS_shapeI_shapeJ_    ;}
    const std::vector<VecDbl>&  GetIntS_shapeI()           const  {return IntS_shapeI_           ;}
    const std::vector<MatVec3>& GetIntS_shapeI_gradshapeJ()const  {return IntS_shapeI_gradshapeJ_;}
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
    std::vector<unsigned int>     quadrature_point_indices_; ///< qp index only
    VecVec3                       qpoints_xyz_             ; ///< qp index only
    std::vector<VecDbl>           shape_value_             ; ///< Node i, then qp
    std::vector<VecVec3>          shape_grad_              ; ///< Node i, then qp
    VecDbl                        JxW_                     ; ///< qp index only
    std::vector<std::vector<int>> face_dof_mappings_       ; ///< Face f,then fi
    size_t                        num_nodes_ = 0;

    bool                          initialized_ = false;

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
    VecVec3                   normals_;                  ///< node i, then qp
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

#endif //CHI_MATH_FINITE_ELEMENT_H