#include "meshcutting.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Raytrace/raytracing.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <algorithm>

//###################################################################
/**Make an edge for a polygon given its edge index.*/
std::pair<uint64_t,uint64_t> chi_mesh::mesh_cutting::
  MakeEdgeFromPolygonEdgeIndex(const chi_mesh::Cell &cell, int edge_index)
{
  int e = edge_index;
  size_t num_verts = cell.vertex_ids.size();

  int next_v = (e < (num_verts-1)) ? e+1 : 0;
  uint64_t v0_id = cell.vertex_ids[e];
  uint64_t v1_id = cell.vertex_ids[next_v];

  return std::make_pair(v0_id,v1_id);
}

//###################################################################
/***/
void chi_mesh::mesh_cutting::
  PopulatePolygonFacesFromVertices(const MeshContinuum &mesh,
                                   const std::vector<uint64_t> &vertex_ids,
                                   chi_mesh::Cell &cell)
{
  cell.faces.clear();
  cell.faces.reserve(vertex_ids.size());
  cell.vertex_ids = vertex_ids;

  size_t num_verts = vertex_ids.size();
  for (size_t v=0; v<num_verts; ++v)
  {
    size_t v1_ref = (v < (num_verts-1))? v+1 : 0;

    uint64_t v0id = cell.vertex_ids[v];
    uint64_t v1id = cell.vertex_ids[v1_ref];

    const chi_mesh::Vertex& v0 = *mesh.vertices[v0id];
    const chi_mesh::Vertex& v1 = *mesh.vertices[v1id];

    chi_mesh::Vector3 v01 = v1 - v0;

    chi_mesh::CellFace face;
    face.vertex_ids = {v0id,v1id};
    face.normal     = v01.Cross(chi_mesh::Normal(0.0,0.0,1.0)).Normalized();
    face.centroid   = (v0+v1)/2.0;

    cell.faces.push_back(face);
  }
}

//###################################################################
/**Cuts a mesh with a plane.*/
void chi_mesh::mesh_cutting::
  CutMeshWithPlane(MeshContinuum& mesh,
                   const Vector3 &plane_point,
                   const Vector3 &plane_normal)
{
  chi_log.Log() << "Cutting mesh with plane.";

  const auto& p = plane_point;
  const auto& n = plane_normal;

  const double merge_tolerance = 1.0e-3;
  const double float_compare   = 1.0e-10;

  /**What the name says.*/
  auto MakeEdgeSetFromPolygonEdgeIndex = [](
    const chi_mesh::Cell& cell,
    int edge_index)
  {
    int e = edge_index;
    size_t num_verts = cell.vertex_ids.size();

    int next_v = (e < (num_verts-1)) ? e+1 : 0;
    uint64_t v0_id = cell.vertex_ids[e];
    uint64_t v1_id = cell.vertex_ids[next_v];

    uint64_t v_max = std::max(v0_id,v1_id);
    uint64_t v_min = std::min(v0_id,v1_id);

    return std::make_pair(v_min,v_max);
  };



  //============================================= Snap vertices to plane within
  //                                              merge tolerance to avoid
  //                                              creating small cells or cutting
  //                                              parallel faces
  // Order N, num_vertices
  size_t num_verts_snapped=0;
  for (auto vertex_ptr : mesh.vertices)
  {
    auto& vertex = *vertex_ptr;
    double d_from_plane = n.Dot(vertex - p);

    if (std::fabs(d_from_plane) < merge_tolerance)
    {
      vertex -= n*d_from_plane;
      ++num_verts_snapped;
    }
  }

  chi_log.Log() << "Number of vertices snapped to plane: "
                << num_verts_snapped;

  //============================================= Determine cells to cut
  // Order N, num_cells
  // A cell is a candidate for cutting if its
  // vertices lay on both sides of the plane.
  // So this algorithm just checks the sense wrt
  // the plane.
  std::vector<chi_mesh::Cell*> cells_to_cut;

  for (auto& cell : mesh.local_cells)
  {
    size_t num_neg_senses = 0;
    size_t num_pos_senses = 0;
    for (auto vid : cell.vertex_ids)
    {
      const auto& x = *mesh.vertices[vid];
      double new_sense = n.Dot(x-p);

      if (new_sense < 0.0) ++num_neg_senses;
      if (new_sense > 0.0) ++num_pos_senses;

      if (num_neg_senses>0 && num_pos_senses>0)
      {
        cells_to_cut.push_back(&cell);
        break;
      }
    }//for vid
  }//for cell
  chi_log.Log() << "Number of cells to cut: " << cells_to_cut.size();

  //============================================= Two-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYGON)
  {
    //====================================== Build unique edges
    size_t num_edges_cut=0;
    std::set<Edge> edges_set;
    for (auto cell_ptr : cells_to_cut)
    {
      const auto& cell = *cell_ptr;
      const size_t num_edges = cell.vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        edges_set.insert(MakeEdgeSetFromPolygonEdgeIndex(cell, e));
    }//for cell - built edges_set

    //====================================== Determine cut edges
    std::vector<ECI> cut_edges;
    {
      for (auto& edge : edges_set)
      {
        const auto& v0 = *mesh.vertices[edge.first];
        const auto& v1 = *mesh.vertices[edge.second];

        chi_mesh::Vector3 cut_point;

        if (CheckPlaneLineIntersect(n,p,v0,v1,cut_point))
        {
          double dv0 = std::fabs((v0-p).Dot(n));
          double dv1 = std::fabs((v1-p).Dot(n));
          if (dv0>float_compare and dv1>float_compare)
          {
            mesh.vertices.push_back(new chi_mesh::Vector3(cut_point));
            cut_edges.emplace_back(edge, mesh.vertices.size() - 1);
            ++num_edges_cut;
          }
        }
      }//for edge - determine cut
    }//populate edges cut

    chi_log.Log() << "Number of cut edges: " << num_edges_cut;

    //====================================== Determine cut vertices
    std::set<uint64_t> cut_vertices;
    {
      for (auto cell_ptr : cells_to_cut)
        for (uint64_t vid : cell_ptr->vertex_ids)
        {
          const auto vertex = *mesh.vertices[vid];
          double dv = std::fabs((vertex-p).Dot(n));
          if (dv<float_compare)
            cut_vertices.insert(vid);
        }//for vid
    }

    chi_log.Log() << "Number of cut vertices: " << cut_vertices.size();

    //====================================== Process cells that are cut
    for (auto cell_ptr : cells_to_cut)
    {
      auto& cell = *(chi_mesh::CellPolygon*)cell_ptr;

      CutPolygon(cut_edges,cut_vertices,p,n,mesh,cell);
    }//for cell_ptr


  }//two-D

  //============================================= Three-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYHEDRON)
  {

  }

  chi_log.Log() << "Done cutting mesh with plane. Num cells = "
                << mesh.local_cells.size();
}