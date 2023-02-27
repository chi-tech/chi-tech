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
    IntV_gradShapeI_gradShapeJ_ = std::move(in_IntV_gradShapeI_gradShapeJ);
    IntV_shapeI_gradshapeJ_     = std::move(in_IntV_shapeI_gradshapeJ    );
    IntV_shapeI_shapeJ_         = std::move(in_IntV_shapeI_shapeJ        );
    IntV_shapeI_                = std::move(in_IntV_shapeI               );
    IntV_gradshapeI_            = std::move(in_IntV_gradshapeI           );
    IntS_shapeI_shapeJ_         = std::move(in_IntS_shapeI_shapeJ    );
    IntS_shapeI_                = std::move(in_IntS_shapeI           );
    IntS_shapeI_gradshapeJ_     = std::move(in_IntS_shapeI_gradshapeJ);
    face_dof_mappings_          = std::move(in_face_dof_mappings);
    num_nodes_                  = in_num_nodes;
  }

  void UnitIntegralData::Reset()
  {
    IntV_gradShapeI_gradShapeJ_.clear();
    IntV_shapeI_gradshapeJ_.clear();
    IntV_shapeI_shapeJ_.clear();
    IntV_shapeI_.clear();
    IntV_gradshapeI_.clear();

    IntS_shapeI_shapeJ_.clear();
    IntS_shapeI_.clear();
    IntS_shapeI_gradshapeJ_.clear();

    face_dof_mappings_.clear();
    num_nodes_=0;
  }

  double UnitIntegralData::
    IntV_gradShapeI_gradShapeJ(unsigned int i,
                               unsigned int j) const
  {
    double value;
    auto& row_I = IntV_gradShapeI_gradShapeJ_.at(i);
    value = row_I.at(j);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    IntV_shapeI_gradshapeJ(unsigned int i,
                           unsigned int j) const
  {
    chi_mesh::Vector3 value;
    auto & row_I = IntV_shapeI_gradshapeJ_.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    IntV_shapeI_shapeJ(unsigned int i,
                       unsigned int j) const
  {
    double value;
    auto& row_I = IntV_shapeI_shapeJ_.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    IntV_shapeI(unsigned int i) const
  {
    double value = IntV_shapeI_.at(i);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    IntV_gradshapeI(unsigned int i) const
  {
    chi_mesh::Vector3 value;
    value = IntV_gradshapeI_.at(i);
    return value;
  }
  double UnitIntegralData::
    IntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j) const
  {
    double value;
    auto& face_data = IntS_shapeI_shapeJ_.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }

  double UnitIntegralData::
    IntS_shapeI(unsigned int face, unsigned int i) const
  {
    double value;
    auto& face_data = IntS_shapeI_.at(face);
    value = face_data.at(i);
    return value;
  }

  chi_mesh::Vector3 UnitIntegralData::
    IntS_shapeI_gradshapeJ(unsigned int face,
                           unsigned int i,
                           unsigned int j) const
  {
    chi_mesh::Vector3 value;
    auto& face_data = IntS_shapeI_gradshapeJ_.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }
}
}

