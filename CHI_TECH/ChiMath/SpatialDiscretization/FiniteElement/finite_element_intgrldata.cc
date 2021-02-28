#include "finite_element.h"

namespace chi_math
{
namespace finite_element
{
  void UnitIntegralData::Initialize(
    MatDbl in_IntV_gradShapeI_gradShapeJ,
    MatVec3 in_IntV_shapeI_gradshapeJ,
    MatDbl in_IntV_shapeI_shapeJ,
    VecDbl in_IntV_shapeI,
    VecVec3 in_IntV_gradshapeI,
    std::vector<MatDbl> in_IntS_shapeI_shapeJ,
    std::vector<VecDbl> in_IntS_shapeI,
    std::vector<MatVec3> in_IntS_shapeI_gradshapeJ,
    std::vector<std::vector<int>> in_face_dof_mappings,
    size_t in_num_nodes)
  {
    m_IntV_gradShapeI_gradShapeJ = std::move(in_IntV_gradShapeI_gradShapeJ);
    m_IntV_shapeI_gradshapeJ     = std::move(in_IntV_shapeI_gradshapeJ    );
    m_IntV_shapeI_shapeJ         = std::move(in_IntV_shapeI_shapeJ        );
    m_IntV_shapeI                = std::move(in_IntV_shapeI               );
    m_IntV_gradshapeI            = std::move(in_IntV_gradshapeI           );
    m_IntS_shapeI_shapeJ         = std::move(in_IntS_shapeI_shapeJ    );
    m_IntS_shapeI                = std::move(in_IntS_shapeI           );
    m_IntS_shapeI_gradshapeJ     = std::move(in_IntS_shapeI_gradshapeJ);
    m_face_dof_mappings          = std::move(in_face_dof_mappings);
    m_num_nodes                  = in_num_nodes;
  }

  void UnitIntegralData::Reset()
  {
    m_IntV_gradShapeI_gradShapeJ.clear();
    m_IntV_shapeI_gradshapeJ.clear();
    m_IntV_shapeI_shapeJ.clear();
    m_IntV_shapeI.clear();
    m_IntV_gradshapeI.clear();

    m_IntS_shapeI_shapeJ.clear();
    m_IntS_shapeI.clear();
    m_IntS_shapeI_gradshapeJ.clear();

    m_face_dof_mappings.clear();
    m_num_nodes=0;
  }

  double UnitIntegralData::
    IntV_gradShapeI_gradShapeJ(unsigned int i,
                               unsigned int j) const
  {
    double value;
    auto& row_I = m_IntV_gradShapeI_gradShapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    IntV_shapeI_gradshapeJ(unsigned int i,
                           unsigned int j) const
  {
    chi_mesh::Vector3 value;
    auto & row_I = m_IntV_shapeI_gradshapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    IntV_shapeI_shapeJ(unsigned int i,
                       unsigned int j) const
  {
    double value;
    auto& row_I = m_IntV_shapeI_shapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    IntV_shapeI(unsigned int i) const
  {
    double value = m_IntV_shapeI.at(i);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    IntV_gradshapeI(unsigned int i) const
  {
    chi_mesh::Vector3 value;
    value = m_IntV_gradshapeI.at(i);
    return value;
  }
  double UnitIntegralData::
    IntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j) const
  {
    double value;
    auto& face_data = m_IntS_shapeI_shapeJ.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }

  double UnitIntegralData::
    IntS_shapeI(unsigned int face, unsigned int i) const
  {
    double value;
    auto& face_data = m_IntS_shapeI.at(face);
    value = face_data.at(i);
    return value;
  }

  chi_mesh::Vector3 UnitIntegralData::
    IntS_shapeI_gradshapeJ(unsigned int face,
                           unsigned int i,
                           unsigned int j) const
  {
    chi_mesh::Vector3 value;
    auto& face_data = m_IntS_shapeI_gradshapeJ.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }
}
}

