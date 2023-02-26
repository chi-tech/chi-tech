#include "cell_mapping_base.h"

#include <utility>

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

chi_math::CellMapping::
  CellMapping(chi_mesh::MeshContinuumConstPtr   in_grid,
              const chi_mesh::Cell& in_cell,
              size_t in_num_nodes,
              std::vector<std::vector<int>> in_face_node_mappings,
              const VandAFunction& volume_area_function) :
  m_grid_ptr(std::move(in_grid)),
  m_cell(in_cell),
  m_num_nodes(in_num_nodes),
  face_node_mappings(std::move(in_face_node_mappings))
{
  volume_area_function(*m_grid_ptr, in_cell, m_volume, m_areas);
}

void chi_math::CellMapping::ComputeCellVolumeAndAreas(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  double& volume,
  std::vector<double>& areas)
{
  switch (cell.Type())
  {
    case chi_mesh::CellType::SLAB:
    {
      const auto& v0 = grid.vertices[cell.vertex_ids_[0]];
      const auto& v1 = grid.vertices[cell.vertex_ids_[1]];

      volume = (v1-v0).Norm();
      areas = {1.0,1.0};
      break;
    }
    case chi_mesh::CellType::POLYGON:
    {
      volume = 0.0;
      const auto& v2 = cell.centroid_;

      size_t num_faces = cell.faces_.size();
      areas.reserve(num_faces);

      for (size_t f=0; f<num_faces; ++f)
      {
        const uint64_t v0i = cell.faces_[f].vertex_ids_[0];
        const uint64_t v1i = cell.faces_[f].vertex_ids_[1];

        const auto& v0 = grid.vertices[v0i];
        const auto& v1 = grid.vertices[v1i];

        areas.push_back((v1-v0).Norm());

        const chi_mesh::Vector3 sidev01 = v1 - v0;
        const chi_mesh::Vector3 sidev02 = v2 - v0;

        double sidedetJ = ((sidev01.x)*(sidev02.y) -
                           (sidev02.x)*(sidev01.y));

        volume += sidedetJ/2.0;
      }//for face

      break;
    }
    case chi_mesh::CellType::POLYHEDRON:
    {
      volume = 0.0;
      const auto& vcc = cell.centroid_;

      size_t num_faces = cell.faces_.size();
      areas.assign(num_faces, 0.0);
      for (size_t f=0; f<num_faces; f++)
      {
        const auto& face = cell.faces_[f];
        const size_t num_edges = face.vertex_ids_.size();
        for (size_t e=0; e<num_edges; ++e)
        {
          size_t ep1 = (e < (num_edges-1))? e+1 : 0;
          uint64_t v0i = face.vertex_ids_[e  ];
          uint64_t v1i = face.vertex_ids_[ep1];

          const auto& v0 = grid.vertices[v0i];
          const auto& v1 = cell.faces_[f].centroid_;
          const auto& v2 = grid.vertices[v1i];
          const auto& v3 = vcc;

          const auto sidev01 = v1-v0;
          const auto sidev02 = v2-v0;
          const auto sidev03 = v3-v0;

          chi_mesh::Matrix3x3 J;

          J.SetColJVec(0,sidev01);
          J.SetColJVec(1,sidev02);
          J.SetColJVec(2,sidev03);

          areas[f] += (sidev01.Cross(sidev02)).Norm()/2.0;
          volume += J.Det()/6.0;
        }//for edge
      }//for face
      break;
    }
    default:
      throw std::logic_error("chi_math::CellMapping::ComputeCellVolume: "
                             "Unsupported cell type.");
  }
}

int chi_math::CellMapping::
  MapFaceNode(size_t face_index, size_t face_node_index) const
{
  try {return face_node_mappings.at(face_index).at(face_node_index);}
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range("chi_math::CellMapping::MapFaceNode: "
          "Either face_index or face_node_index is out of range");
  }
}

void
chi_math::CellMapping::
  InitializeAllQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data,
  std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data) const
{
  InitializeVolumeQuadraturePointData(internal_data);
  faces_qp_data.resize(face_node_mappings.size());
  for (size_t f = 0; f < faces_qp_data.size(); ++f)
    InitializeFaceQuadraturePointData(f, faces_qp_data[f]);
}


void
chi_math::CellMapping::
  ComputeUnitIntegrals(
    chi_math::finite_element::UnitIntegralData& ui_data) const
{
  //  quadrature point data
  chi_math::finite_element::InternalQuadraturePointData internal_data;
  std::vector<chi_math::finite_element::FaceQuadraturePointData> faces_qp_data;
  InitializeAllQuadraturePointData(internal_data, faces_qp_data);

  //  integrals
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
                     face_node_mappings,
                     m_num_nodes);
}


chi_math::finite_element::InternalQuadraturePointData
  chi_math::CellMapping::MakeVolumeQuadraturePointData() const
{
  chi_math::finite_element::InternalQuadraturePointData qp_data;
  InitializeVolumeQuadraturePointData(qp_data);
  return qp_data;
}

chi_math::finite_element::FaceQuadraturePointData
chi_math::CellMapping::MakeFaceQuadraturePointData(size_t face_index) const
{
  chi_math::finite_element::FaceQuadraturePointData qp_data;
  InitializeFaceQuadraturePointData(face_index, qp_data);
  return qp_data;
}