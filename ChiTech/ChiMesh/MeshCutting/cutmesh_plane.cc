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
  const double float_compare   = 1.0e-18;

  /**Scans a vector of EdgeCutInfo to see if it contains the given edge.*/
  auto VecOfECIHasEdge = [](
    const std::vector<ECI>& vec_edge_cut_info,
    const Edge& edge)
  {
    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto find_result =  std::find_if(vec_edge_cut_info.begin(),
                                     vec_edge_cut_info.end(),
                                     std::bind(Comparator,Arg1,edge));

    if (find_result != vec_edge_cut_info.end())
      return std::make_pair(true,*find_result);

    return std::make_pair(false,*find_result);
  };

  /**Scans a vector of EdgeCutInfo to see if it contains the given edge.*/
  auto VecOfECIGetInfoMatchingEdge = [](
    const std::vector<ECI>& vec_edge_cut_info,
    const Edge& edge)
  {
    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto find_result =  std::find_if(vec_edge_cut_info.begin(),
                                     vec_edge_cut_info.end(),
                                     std::bind(Comparator,Arg1,edge));

    if (find_result != vec_edge_cut_info.end())
      return *find_result;
    else
      throw std::logic_error("CutMeshWithPlane:VecECIGetInfoMatchingEdge");
  };

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

  /**What the name says.*/
  auto MakeEdgeFromPolygonEdgeIndex = [](
    const chi_mesh::Cell& cell,
    int edge_index)
  {
    int e = edge_index;
    size_t num_verts = cell.vertex_ids.size();

    int next_v = (e < (num_verts-1)) ? e+1 : 0;
    uint64_t v0_id = cell.vertex_ids[e];
    uint64_t v1_id = cell.vertex_ids[next_v];

    return std::make_pair(v0_id,v1_id);
  };

  /**Makes a complete polygon from just the vertex ids.*/
  auto PopulatePolygonFacesFromVertices = [&mesh](
    const std::vector<uint64_t>& vertex_ids,
    chi_mesh::Cell& cell)
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

    chi_log.Log() << "Number of cut vertices: " << num_edges_cut;

    //====================================== Process cells that are cut
    for (auto cell_ptr : cells_to_cut)
    {
      auto& cell = *(chi_mesh::CellPolygon*)cell_ptr;
      const size_t num_edges = cell.vertex_ids.size();

      //=============================== Find the first edge of this cell
      //                                that is a cut-edge
      // If no edge is found then this
      // could indicate a cut through vertices
      int first_cut_edge_index =-1;
      ECI first_cut_ECI_info;
      for (int e=0; e<num_edges; ++e)
      {
        auto edge = MakeEdgeSetFromPolygonEdgeIndex(cell, e);
        auto result = VecOfECIHasEdge(cut_edges, edge);
        if (result.first)
        {
          first_cut_edge_index = e;
          first_cut_ECI_info = result.second;
          break;
        }
      }//for edge e

      //=============================== Process case where no edge is cut
      if (first_cut_edge_index<0)
      {
        std::vector<uint64_t> cur_cell_vertex_ids;
        std::vector<uint64_t> new_cell_vertex_ids;
        bool on_current_cell = true;
        for (uint64_t vid : cell.vertex_ids)
        {
          const auto& vertex = *mesh.vertices[vid];
          if (std::fabs((vertex-p).Dot(n)) < float_compare)
          {
            if (on_current_cell)
            {
              cur_cell_vertex_ids.push_back(vid);
              on_current_cell = false;
            }
            else
            {
              new_cell_vertex_ids.push_back(vid);
              on_current_cell = true;
            }
          }

          if (on_current_cell) cur_cell_vertex_ids.push_back(vid);
          else                 new_cell_vertex_ids.push_back(vid);
        }//for vid

        auto new_cell = new chi_mesh::CellPolygon;

        PopulatePolygonFacesFromVertices(cur_cell_vertex_ids,cell);
        PopulatePolygonFacesFromVertices(new_cell_vertex_ids,*new_cell);
        mesh.cells.push_back(new_cell);
      }//if no-edges cut
      else
      {
        //=============================== Find the second edge this cell
        //                                that is counter clock wise from
        //                                the first one found
        int second_cut_edge_index = -1;
        ECI second_cut_ECI_info;
        for (int eref=0,e=first_cut_edge_index; eref<(num_edges-1); ++eref)
        {
          e = (e<(num_edges-1))? e+1 : 0;

          auto edge = MakeEdgeSetFromPolygonEdgeIndex(cell, e);
          auto result = VecOfECIHasEdge(cut_edges, edge);
          if (result.first)
          {
            second_cut_edge_index = e;
            second_cut_ECI_info = result.second;
            break;
          }
        }

        if (second_cut_edge_index < 0) throw std::logic_error("second_cut_edge<0");

        //=============================== Making cur_cell_vertex id changes
        std::vector<uint64_t> cur_cell_vertex_ids;
        {
          int e = second_cut_edge_index;
          while (e != first_cut_edge_index)
          {
            Edge edge_set = MakeEdgeSetFromPolygonEdgeIndex(cell, e);

            if (e == second_cut_edge_index)
            {
              ECI cut_info = VecOfECIGetInfoMatchingEdge(cut_edges, edge_set);
              cur_cell_vertex_ids.push_back(cut_info.cut_point_id);
            }
            else
            {
              Edge edge = MakeEdgeFromPolygonEdgeIndex(cell,e);
              cur_cell_vertex_ids.push_back(edge.first);
            }

            e = (e<(num_edges-1))? e+1 : 0;
          }

          Edge edge = MakeEdgeFromPolygonEdgeIndex(cell,e);
          cur_cell_vertex_ids.push_back(edge.first);
          cur_cell_vertex_ids.push_back(first_cut_ECI_info.cut_point_id);
        }

        //=============================== Making new cell vertex ids
        std::vector<uint64_t> new_cell_vertex_ids;
        {
          int e = first_cut_edge_index;
          while (e != second_cut_edge_index)
          {
            Edge edge_set = MakeEdgeSetFromPolygonEdgeIndex(cell, e);

            if (e == first_cut_edge_index)
            {
              ECI cut_info = VecOfECIGetInfoMatchingEdge(cut_edges, edge_set);
              new_cell_vertex_ids.push_back(cut_info.cut_point_id);
            }
            else
            {
              Edge edge = MakeEdgeFromPolygonEdgeIndex(cell,e);
              new_cell_vertex_ids.push_back(edge.first);
            }

            e = (e<(num_edges-1))? e+1 : 0;
          }

          Edge edge = MakeEdgeFromPolygonEdgeIndex(cell,e);
          new_cell_vertex_ids.push_back(edge.first);
          new_cell_vertex_ids.push_back(second_cut_ECI_info.cut_point_id);
        }

        auto new_cell = new chi_mesh::CellPolygon;

        PopulatePolygonFacesFromVertices(cur_cell_vertex_ids,cell);
        PopulatePolygonFacesFromVertices(new_cell_vertex_ids,*new_cell);
        mesh.cells.push_back(new_cell);

        chi_log.Log() << "current cell:";
        for (auto vid : cell.vertex_ids)
          chi_log.Log() << vid;

        chi_log.Log() << "new cell:";
        for (auto vid : new_cell->vertex_ids)
          chi_log.Log() << vid;
      }
    }//for cell_ptr


  }//two-D

  //============================================= Three-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYHEDRON)
  {

  }

  chi_log.Log() << "Done cutting mesh with plane. Num cells = "
                << mesh.local_cells.size();
}