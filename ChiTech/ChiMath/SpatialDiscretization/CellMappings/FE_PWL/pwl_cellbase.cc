#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_log.h"

std::vector<chi_mesh::Vector3> chi_math::CellMappingFE_PWL::GetNodeLocations() const
{
  return m_node_locations;
}

/** This section just determines a mapping of face dofs
to cell dofs. This is pretty simple since we can
just loop over each face dof then subsequently
loop over cell dofs, if the face dof node index equals
the cell dof node index then the mapping is assigned.

This mapping is not used by any of the methods in
this class but is used by methods requiring the
surface integrals of the shape functions.*/
std::vector<std::vector<int>> chi_math::CellMappingFE_PWL::
  MakeFaceNodeMapping(const chi_mesh::Cell &cell)
{
  const size_t num_faces = cell.faces.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (auto& face : cell.faces)
  {
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face.vertex_ids.size());
    for (uint64_t fvid : face.vertex_ids)
    {
      int mapping = -1;
      for (size_t ci=0; ci<cell.vertex_ids.size(); ci++)
      {
        if (fvid == cell.vertex_ids[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      }//for cell i
      if (mapping<0)
      {
        chi::log.LogAllError() << "Unknown face mapping encountered. "
                                  "pwl_polyhedron.h";
        chi::Exit(EXIT_FAILURE);
      }
      face_dof_mapping.push_back(mapping);
    }//for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

void
chi_math::CellMappingFE_PWL::
  ComputeWeightedUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const
{
  //  quadrature point data
  chi_math::finite_element::InternalQuadraturePointData internal_data;
  std::vector<chi_math::finite_element::FaceQuadraturePointData> faces_qp_data;
  InitializeAllQuadraturePointData(internal_data, faces_qp_data);


  //  weighted JxW
  std::vector<double> V_JxW;
  for (const auto& qp : internal_data.QuadraturePointIndices())
  {
    const auto jxw =
      SpatialWeightFunction(internal_data.QPointXYZ(qp)) *
      internal_data.JxW(qp);
    V_JxW.emplace_back(jxw);
  }

  std::vector<std::vector<double>> F_JxW(faces_qp_data.size());
  for (size_t f = 0; f < faces_qp_data.size(); ++f)
    for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
    {
      const auto jxw =
        SpatialWeightFunction(faces_qp_data[f].QPointXYZ(qp)) *
        faces_qp_data[f].JxW(qp);
      F_JxW[f].emplace_back(jxw);
    }


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
             V_JxW[qp];

        IntV_shapeI_gradshapeJ[i][j]
          += internal_data.ShapeValue(i, qp) *
             internal_data.ShapeGrad(j, qp) *
             V_JxW[qp];

        IntV_shapeI_shapeJ[i][j]
          += internal_data.ShapeValue(i, qp) *
             internal_data.ShapeValue(j, qp) *
             V_JxW[qp];
      }// for qp
    }// for j

    for (const auto& qp : internal_data.QuadraturePointIndices())
    {
      IntV_shapeI[i]
        += internal_data.ShapeValue(i, qp) * V_JxW[qp];

      IntV_gradshapeI[i]
        += internal_data.ShapeGrad(i, qp) * V_JxW[qp];
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
               F_JxW[f][qp];

          IntS_shapeI_gradshapeJ[f][i][j]
            += faces_qp_data[f].ShapeValue(i, qp) *
               faces_qp_data[f].ShapeGrad(j, qp) *
               F_JxW[f][qp];
        }// for qp
      }//for j

      for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
      {
        IntS_shapeI[f][i]
          += faces_qp_data[f].ShapeValue(i, qp) * F_JxW[f][qp];
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
                     face_node_mappings,
                     m_num_nodes);
}

std::vector<chi_mesh::Vector3>
chi_math::CellMappingFE_PWL::
  GetVertexLocations(const chi_mesh::MeshContinuum &grid,
                     const chi_mesh::Cell &cell)
{
  std::vector<chi_mesh::Vector3> verts;
  verts.reserve(cell.vertex_ids.size());

  for (const auto vid : cell.vertex_ids)
    verts.push_back(grid.vertices[vid]);

  return verts;
}