#include "meshcutting.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <algorithm>

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

  typedef std::pair<uint64_t, uint64_t> Edge;
  typedef chi_mesh::Vector3 Vec3;

  struct EdgeCutInfo
  {
    Edge vertex_ids;
    uint64_t cut_point_id;

    explicit EdgeCutInfo(Edge in_edge,
                         uint64_t in_cutpoint_id) :
      vertex_ids(std::move(in_edge)),
      cut_point_id(in_cutpoint_id)
    {}

    static
    bool Comparator(const EdgeCutInfo& edge_cut_info,
                    const Edge& ref_edge)
    {
      return ref_edge == edge_cut_info.vertex_ids;
    }
  };

  auto VecEdgeCutInfoHasEdge = [](const std::vector<EdgeCutInfo>& vec_edge_cut_info,
                                  const Edge& edge)
  {
    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &EdgeCutInfo::Comparator;

    auto find_result =  std::find_if(vec_edge_cut_info.begin(),
                                     vec_edge_cut_info.end(),
                                     std::bind(Comparator,Arg1,edge));

    if (find_result != vec_edge_cut_info.end())
      return true;

    return false;
  };

  auto MakeEdgeFromPolygonEdgeIndex = [](const chi_mesh::Cell& cell,
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
  for (auto vertex_ptr : mesh.vertices)
  {
    auto& vertex = *vertex_ptr;
    double d_from_plane = n.Dot(vertex - p);

    if (std::fabs(d_from_plane) < merge_tolerance)
      vertex -= n*d_from_plane;
  }

  //============================================= Determine cells to cut
  // Order N, num_cells
  std::vector<chi_mesh::Cell*> cells_to_cut;

  for (auto& cell : mesh.local_cells)
  {
    double old_sense = 0.0;
    for (auto vid : cell.vertex_ids)
    {
      const auto& x = *mesh.vertices[vid];
      double new_sense = n.Dot(x-p);

      if (new_sense*old_sense < 0.0)
      {
        cells_to_cut.push_back(&cell);
        break;
      }
      old_sense = new_sense;
    }//for vid
  }//for cell
  chi_log.Log() << "Num cells to cut: " << cells_to_cut.size();

  // Facts:
  // - The fundamental unit that can be cut by a plane
  //   is a straight edge.
  // - An edge may be shared by multiple cells.
  // - Cutting through a cell can create a cut-face or a
  //   a face that is completely removed from a cell.
  // - For a polygon an edge is also a face.
  // - For a polyhedron an edge may be shared by multiple faces.

  //============================================= Two-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYGON)
  {
    //====================================== Build unique edges
    size_t num_edges_cut=0;
    std::set<Edge> edges_set;
    for (auto cell_ptr : cells_to_cut)
    {
      auto& cell = *cell_ptr;
      const size_t num_edges = cell.vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        edges_set.insert(MakeEdgeFromPolygonEdgeIndex(cell,e));
    }//for cell - built edges_set

    chi_log.Log() << "Edges to be inspected for cutting:";
    for (const auto& edge : edges_set)
      chi_log.Log() << edge.first << " " << edge.second;

    //====================================== Determine edges to be cut
    std::vector<EdgeCutInfo> edges_to_cut;
    {
      for (auto& edge : edges_set)
      {
        const auto& v0 = *mesh.vertices[edge.first];
        const auto& v1 = *mesh.vertices[edge.second];

        chi_mesh::Vector3 cut_point;

        if (CheckPlaneLineIntersect(n,p,v0,v1,cut_point))
        {
          mesh.vertices.push_back(new chi_mesh::Vector3(cut_point));
          edges_to_cut.emplace_back(edge,mesh.vertices.size()-1);
          ++num_edges_cut;
        }
      }//for edge - determine cut
    }//populate edges cut

    chi_log.Log() << "Num edges to cut: " << num_edges_cut;

    //====================================== Process cells that are cut
    for (auto cell_ptr : cells_to_cut)
    {
      auto& cell = *cell_ptr;
      const size_t num_edges = cell.vertex_ids.size();

      //=============================== Find the first edge of this cell
      //                                that is a cut-edge
      // If no edge is found then this is
      // bad and the code must puke.
      int      first_cut_edge_index =-1;
      for (int e=0; e<num_edges; ++e)
      {
        auto edge = MakeEdgeFromPolygonEdgeIndex(cell,e);

        if (VecEdgeCutInfoHasEdge(edges_to_cut, edge))
        {
          first_cut_edge_index = e;
          break;
        }
      }//for edge e
      if (first_cut_edge_index < 0)
        throw std::logic_error(std::string(__FUNCTION__) + " line " +
                               std::to_string(__LINE__));

      //=============================== Trace from this edge, ccw, to another
      //                                edge that is also cut
      // If no second cut is found then
      // it is ok. Only one edge is split.
      int second_cut_edge_index=-1;
      {
        int e = (first_cut_edge_index<(num_edges-1))? first_cut_edge_index+1 : 0;
        while (e != first_cut_edge_index)
        {
          auto edge = MakeEdgeFromPolygonEdgeIndex(cell,e);

          if (VecEdgeCutInfoHasEdge(edges_to_cut, edge))
          {
            second_cut_edge_index = e;
            break;
          }

          e = (e<(num_edges-1))? e+1 : 0;
        }
      }//second edge logic

      //=============================== Second not found
      if (second_cut_edge_index < 0)
      {
        std::vector<CellFace> new_faces;
        for (int e=0; e<num_edges; ++e)
        {
          if (e == first_cut_edge_index)
          {
            CellFace new_face;
          }
        }

      }//second edge cell modification


    }//for cell_ptr


  }//two-D

  //============================================= Three-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYHEDRON)
  {

  }

  chi_log.Log() << "Done cutting mesh with plane.";
}