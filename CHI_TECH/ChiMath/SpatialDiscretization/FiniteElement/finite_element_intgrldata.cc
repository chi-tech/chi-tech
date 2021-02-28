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
    IntV_gradShapeI_gradShapeJ = std::move(in_IntV_gradShapeI_gradShapeJ);
    IntV_shapeI_gradshapeJ     = std::move(in_IntV_shapeI_gradshapeJ    );
    IntV_shapeI_shapeJ         = std::move(in_IntV_shapeI_shapeJ        );
    IntV_shapeI                = std::move(in_IntV_shapeI               );
    IntV_gradshapeI            = std::move(in_IntV_gradshapeI           );
    IntS_shapeI_shapeJ         = std::move(in_IntS_shapeI_shapeJ    );
    IntS_shapeI                = std::move(in_IntS_shapeI           );
    IntS_shapeI_gradshapeJ     = std::move(in_IntS_shapeI_gradshapeJ);
    face_dof_mappings          = std::move(in_face_dof_mappings);
    num_nodes                  = in_num_nodes;
  }

  void UnitIntegralData::Reset()
  {
    IntV_gradShapeI_gradShapeJ.clear();
    IntV_shapeI_gradshapeJ.clear();
    IntV_shapeI_shapeJ.clear();
    IntV_shapeI.clear();
    IntV_gradshapeI.clear();

    IntS_shapeI_shapeJ.clear();
    IntS_shapeI.clear();
    IntS_shapeI_gradshapeJ.clear();

    face_dof_mappings.clear();
    num_nodes=0;
  }

  double UnitIntegralData::
    FIntV_gradShapeI_gradShapeJ(unsigned int i,
                                unsigned int j)
  {
    double value;
    auto& row_I = IntV_gradShapeI_gradShapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    FIntV_shapeI_gradshapeJ(unsigned int i,
                            unsigned int j)
  {
    chi_mesh::Vector3 value;
    auto & row_I = IntV_shapeI_gradshapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    FIntV_shapeI_shapeJ(unsigned int i,
                        unsigned int j)
  {
    double value;
    auto& row_I = IntV_shapeI_shapeJ.at(i);
    value = row_I.at(j);
    return value;
  }
  double UnitIntegralData::
    FIntV_shapeI(unsigned int i)
  {
    double value = IntV_shapeI.at(i);
    return value;
  }
  chi_mesh::Vector3 UnitIntegralData::
    FIntV_gradshapeI(unsigned int i)
  {
    chi_mesh::Vector3 value;
    value = IntV_gradshapeI.at(i);
    return value;
  }
  double UnitIntegralData::
    FIntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j)
  {
    double value;
    auto& face_data = IntS_shapeI_shapeJ.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }

  double UnitIntegralData::
    FIntS_shapeI(unsigned int face, unsigned int i)
  {
    double value;
    auto& face_data = IntS_shapeI.at(i);
    value = face_data.at(face);
    return value;
  }

  chi_mesh::Vector3 UnitIntegralData::
    FIntS_shapeI_gradshapeJ(unsigned int face,
                            unsigned int i,
                            unsigned int j)
  {
    chi_mesh::Vector3 value;
    auto& face_data = IntS_shapeI_gradshapeJ.at(face);
    auto& rowI      = face_data.at(i);
    value = rowI.at(j);
    return value;
  }
}
}

