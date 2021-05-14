#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"


void
CellMappingFE_PWL::
  InitializeAllQuadraturePointData(
    chi_math::finite_element::InternalQuadraturePointData& internal_data,
    std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) const
{
  InitializeVolumeQuadraturePointData(internal_data);
  faces_qp_data.resize(face_dof_mappings.size());
  for (size_t f = 0; f < faces_qp_data.size(); ++f)
    InitializeFaceQuadraturePointData(f, faces_qp_data[f]);
}


void
CellMappingFE_PWL::
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const
{
  //  quadrature point data
  chi_math::finite_element::InternalQuadraturePointData internal_data;
  std::vector<chi_math::finite_element::FaceQuadraturePointData> faces_qp_data;
  InitializeAllQuadraturePointData(internal_data, faces_qp_data);


  //  integrals
  using VecDbl  = std::vector<double>;
  using MatDbl  = std::vector<VecDbl>;
  using VecVec3 = std::vector<chi_mesh::Vector3>;
  using MatVec3 = std::vector<VecVec3>;

  const auto n_dof_per_cell = internal_data.NumNodes();

  MatDbl  IntV_gradshapeI_gradshapeJ(n_dof_per_cell, VecDbl(n_dof_per_cell));
  MatVec3 IntV_shapeI_gradshapeJ(n_dof_per_cell, VecVec3(n_dof_per_cell));
  MatDbl  IntV_shapeI_shapeJ(n_dof_per_cell, VecDbl(n_dof_per_cell));
  VecDbl  IntV_shapeI(n_dof_per_cell);
  VecVec3 IntV_gradshapeI(n_dof_per_cell);

  std::vector<MatDbl>  IntS_shapeI_shapeJ(faces_qp_data.size());
  std::vector<VecDbl>  IntS_shapeI(faces_qp_data.size());
  std::vector<MatVec3> IntS_shapeI_gradshapeJ(faces_qp_data.size());


  //  volume integrals
  for (unsigned int i = 0; i < n_dof_per_cell; ++i)
  {
    for (unsigned int j = 0; j < n_dof_per_cell; ++j)
    {
      for (const auto& qp : internal_data.QuadraturePointIndices())
      {
        IntV_gradshapeI_gradshapeJ[i][j]
          += internal_data.ShapeGrad(i, qp).Dot(internal_data.ShapeGrad(j, qp)) *
             internal_data.JxW(qp);

        IntV_shapeI_gradshapeJ[i][j]
          += internal_data.ShapeValue(i, qp) *
             internal_data.ShapeGrad(j, qp) *
             internal_data.JxW(qp);

        IntV_shapeI_shapeJ[i][j]
          += internal_data.ShapeValue(i, qp) *
             internal_data.ShapeValue(j, qp) *
             internal_data.JxW(qp);
      }// for qp
    }// for j

    for (const auto& qp : internal_data.QuadraturePointIndices())
    {
      IntV_shapeI[i]
        += internal_data.ShapeValue(i, qp) * internal_data.JxW(qp);

      IntV_gradshapeI[i]
        += internal_data.ShapeGrad(i, qp) * internal_data.JxW(qp);
    }// for qp
  }//for i


  //  surface integrals
  for (size_t f = 0; f < faces_qp_data.size(); ++f)
  {
    IntS_shapeI_shapeJ[f].resize(n_dof_per_cell, VecDbl(n_dof_per_cell));
    IntS_shapeI[f].resize(n_dof_per_cell);
    IntS_shapeI_gradshapeJ[f].resize(n_dof_per_cell, VecVec3(n_dof_per_cell));

    for (unsigned int i = 0; i < n_dof_per_cell; ++i)
    {
      for (unsigned int j = 0; j < n_dof_per_cell; ++j)
      {
        for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
        {
          IntS_shapeI_shapeJ[f][i][j]
            += faces_qp_data[f].ShapeValue(i, qp) *
               faces_qp_data[f].ShapeValue(j, qp) *
               faces_qp_data[f].JxW(qp);

          IntS_shapeI_gradshapeJ[f][i][j]
            += faces_qp_data[f].ShapeValue(i, qp) *
               faces_qp_data[f].ShapeGrad(j, qp) *
               faces_qp_data[f].JxW(qp);
        }// for qp
      }//for j

      for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
      {
        IntS_shapeI[f][i]
          += faces_qp_data[f].ShapeValue(i, qp) * faces_qp_data[f].JxW(qp);
      }// for qp
    }//for i
  }//for f


  //  unit integral data
  ui_data.Initialize(IntV_gradshapeI_gradshapeJ,
                     IntV_shapeI_gradshapeJ,
                     IntV_shapeI_shapeJ,
                     IntV_shapeI,
                     IntV_gradshapeI,
                     IntS_shapeI_shapeJ,
                     IntS_shapeI,
                     IntS_shapeI_gradshapeJ,
                     face_dof_mappings,
                     num_nodes);
}
