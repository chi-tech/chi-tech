#include "meshcutting.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <algorithm>
#include <unistd.h>

//###################################################################
/**Cuts a polygon.*/
void chi_mesh::mesh_cutting::
  CutPolygon(const std::vector<ECI> &cut_edges,
             const std::set<uint64_t> &cut_vertices,
             const Vector3 &plane_point,
             const Vector3 &plane_normal,
             MeshContinuum &mesh,
             chi_mesh::CellPolygon &cell)
{
//  chi_log.Log() << "Processing cell cut";

  const auto& p = plane_point;
  const auto& n = plane_normal;

  const double merge_tolerance = 1.0e-3;
  const double float_compare   = 1.0e-10;

  auto VertexIsCut = [&cut_vertices](uint64_t vid)
  {
    auto result = cut_vertices.find(vid);

    if (result != cut_vertices.end())
      return true;

    return false;
  };

  auto EdgeIsCut = [&cut_edges](const Edge& edge)
  {
    Edge edge_set(std::min(edge.first,edge.second),
                  std::max(edge.first,edge.second));

    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto result = std::find_if(cut_edges.begin(),
                               cut_edges.end(),
                               std::bind(Comparator,Arg1,edge_set));

    if (result != cut_edges.end())
      return std::make_pair(true,*result);

    return std::make_pair(false,*result);
  };

  //============================================= Set vertex and edge cut flags
  size_t num_verts = cell.vertex_ids.size();
  size_t num_edges = num_verts;

  std::vector<bool> vertex_cut_flags(num_verts,false);
  std::vector<bool> edge_cut_flags(num_edges,false);
  std::vector<ECI>  edge_cut_info(num_edges);

  for (size_t e=0; e<num_edges; ++e)
  {
    vertex_cut_flags[e] = VertexIsCut(cell.vertex_ids[e]);

    auto edge = MakeEdgeFromPolygonEdgeIndex(cell, e);
    auto cut_nature = EdgeIsCut(edge);
    edge_cut_flags[e] = cut_nature.first;
    if (cut_nature.first) edge_cut_info[e] = cut_nature.second;

    edge_cut_info[e].vertex_ids = edge;
  }

  //============================================= Lamda for edge loop
  enum class CurVertex
  {
    AT_FIRST,
    AT_CUT_POINT,
    AT_SECOND,
    NONE
  };

  struct CurCutInfo
  {
    int which_edge=0;
    CurVertex which_vertex=CurVertex::AT_FIRST;

    CurCutInfo() = default;

    CurCutInfo(int in_which_edge, CurVertex in_which_vertex) :
      which_edge(in_which_edge),
      which_vertex(in_which_vertex) {}
  };

  /**Follow edge loop.*/
  auto GetVerticesTillNextCut =
    [&cell,&edge_cut_flags,&edge_cut_info,&VertexIsCut](
      CurCutInfo start_cut_info)
  {
    size_t num_verts = cell.vertex_ids.size();
    std::vector<uint64_t> vertex_ids;
    vertex_ids.reserve(num_verts);

    int e = start_cut_info.which_edge;

    auto end_type = CurVertex::NONE;

    switch (start_cut_info.which_vertex)
    {
      case CurVertex::AT_FIRST:
      {
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.first);
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);

        if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
        {
          end_type = CurVertex::AT_SECOND;
          goto skip_to_return_portion;
        }

        break;
      }
      case CurVertex::AT_CUT_POINT:
      {
        vertex_ids.push_back(edge_cut_info[e].cut_point_id);
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);

        if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
        {
          end_type = CurVertex::AT_SECOND;
          goto skip_to_return_portion;
        }

        break;
      }
      case CurVertex::AT_SECOND:
      case CurVertex::NONE:
      default:
        break;
    }

    //Look at downstream ccw edges and check for
    //edges cut or end-point cuts
    for (int eref=0; eref<num_verts; ++eref)
    {
      e = (e<(num_verts-1))? e+1 : 0;

      if (e == start_cut_info.which_edge)
      {
        break;
      }

      if (edge_cut_flags[e])
      {
        vertex_ids.push_back(edge_cut_info[e].cut_point_id);
        end_type = CurVertex::AT_CUT_POINT;
        break;
      }
      else if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
      {
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);
        end_type = CurVertex::AT_SECOND;
        break;
      }
      else
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);
    }

  skip_to_return_portion:
    CurCutInfo end_cut_info(e, end_type);

    return std::make_pair(vertex_ids,end_cut_info);
  };

  typedef std::pair<std::vector<uint64_t>,CurCutInfo> LoopInfo;

  //============================================= Process all edges
  std::vector<chi_mesh::Cell*> cells_to_add_to_mesh;
  std::vector<CurCutInfo> cut_history;
  for (size_t e=0; e<num_edges; ++e)
  {
//    chi_log.Log() << "Inspecting edge " << e
//                  << " from " << edge_cut_info[e].vertex_ids.first
//                  << " to " << edge_cut_info[e].vertex_ids.second;

    std::vector<uint64_t> verts_to_next_cut;
    int           end_e = e;
    CurVertex end_type = CurVertex::NONE;
    LoopInfo loop_info;

    if (vertex_cut_flags[e])
    {
//      chi_log.Log() << "  Vertex cut " << e;
      cut_history.emplace_back(e,CurVertex::AT_FIRST);
      loop_info = GetVerticesTillNextCut(cut_history.back());
    }
    else if (edge_cut_flags[e])
    {
//      chi_log.Log() << "  Edge cut " << e;
      cut_history.emplace_back(e,CurVertex::AT_CUT_POINT);
      loop_info = GetVerticesTillNextCut(cut_history.back());
    }
    else continue;

    verts_to_next_cut = loop_info.first;
    end_e = loop_info.second.which_edge;
    end_type = loop_info.second.which_vertex;

//    chi_log.Log() << "  edge_end & num_verts " << end_e
//                  << " " << verts_to_next_cut.size();

    if (end_e < e) //if looped past num_edges
      e = (int)num_edges-1; //stop search
    else if (end_type == CurVertex::AT_SECOND)
      e = end_e;            //resume search after end_e
    else
      e = end_e - 1;        //resume search at end_e

    auto new_cell = new chi_mesh::CellPolygon;
    PopulatePolygonFacesFromVertices(mesh,verts_to_next_cut,*new_cell);
    cells_to_add_to_mesh.push_back(new_cell);
  }//for e

//  chi_log.Log() << "Cells to add to the mesh: " << cells_to_add_to_mesh.size();

  if (cells_to_add_to_mesh.size()>2)
  {
//    //=========================================== Make cell-edge pairs, with the
//    //                                            edge being on the cut-plane
//    typedef std::pair<chi_mesh::Cell*,int> CellEdgePair;
//    std::vector<CellEdgePair> cell_edge_pairs_on_plane;
//    size_t cell_counter = 0;
//    for (auto cell_ptr : cells_to_add_to_mesh)
//    {
//      for (int e=0; e<cell_ptr->vertex_ids.size(); ++e)
//      {
//        auto edge = MakeEdgeFromPolygonEdgeIndex(*cell_ptr,e);
//
//        const auto& v0 = *mesh.vertices[edge.first];
//        const auto& v1 = *mesh.vertices[edge.second];
//
//        double dv0 = (v0-p).Dot(n);
//        double dv1 = (v1-p).Dot(n);
//
//        bool v0_on_plane = std::fabs(dv0) < float_compare;
//        bool v1_on_plane = std::fabs(dv1) < float_compare;
//
//        chi_log.Log() << cell_counter << " "
//                      << edge.first << "->" << edge.second << " "
//                      << dv0 << " " << dv1;
//
//        if (v0_on_plane and v1_on_plane)
//          cell_edge_pairs_on_plane.emplace_back(cell_ptr,e);
//      }
//      ++cell_counter;
//    }
//
//    chi_log.Log() << "Number of edges on plane: " << cell_edge_pairs_on_plane.size();
//
//    //=========================================== Find an edge of a cell-edge
//    //                                            that contains the other edges
//    CellEdgePair master_cell_edge_pair(nullptr,-1);
//    for (const auto& cell_edge_i_pair : cell_edge_pairs_on_plane)
//    {
//      bool assumption_contains_all = true;
//      const auto& cell_i = *cell_edge_i_pair.first;
//      auto edge_i = MakeEdgeFromPolygonEdgeIndex(cell_i, cell_edge_i_pair.second);
//
//      const auto& v0_edge_i = *mesh.vertices[edge_i.first];
//      const auto& v1_edge_i = *mesh.vertices[edge_i.second];
//
//      auto   v01_edge_i = v1_edge_i - v0_edge_i;
//      double D = v01_edge_i.Norm();
//      auto   hat_v01_edge_i = v01_edge_i / D;
//
//      for (const auto& cell_edge_j_pair : cell_edge_pairs_on_plane)
//        if (cell_edge_i_pair.first == cell_edge_j_pair.first)
//          continue;
//        else
//        {
//          const auto& cell_j = *cell_edge_j_pair.first;
//          auto edge_j = MakeEdgeFromPolygonEdgeIndex(cell_j, cell_edge_j_pair.second);
//
//          const auto& v0_edge_j = *mesh.vertices[edge_j.first];
//          const auto& v1_edge_j = *mesh.vertices[edge_j.second];
//
//          auto v0j_v0i = v0_edge_j - v0_edge_i;
//          auto v1j_v0i = v1_edge_j - v0_edge_i;
//
//          double d_v0j_v0i = v0j_v0i.Dot(hat_v01_edge_i);
//          double d_v1j_v0i = v1j_v0i.Dot(hat_v01_edge_i);
//
//          bool contains_v0j = (d_v0j_v0i >= (0.0-float_compare)) and
//                              (d_v0j_v0i < (D+float_compare));
//          bool contains_v1j = (d_v1j_v0i >= (0.0-float_compare)) and
//                              (d_v1j_v0i < (D+float_compare));
//
//          if ((not contains_v0j) or (not contains_v1j))
//          {
//            assumption_contains_all = false;
//            break;
//          }
//        }
//
//      if (assumption_contains_all)
//      {
//        master_cell_edge_pair = cell_edge_i_pair;
//        break;
//      }
//    }//for cell-edge pair
//
//    if (master_cell_edge_pair.first != nullptr)
//    {
//      std::set<uint64_t> vertex_ids_set;
//
//      auto& master_cell = *master_cell_edge_pair.first;
//      int   master_edge_index = master_cell_edge_pair.second;
//      auto master_edge = MakeEdgeFromPolygonEdgeIndex(master_cell,master_edge_index);
//      int master_vid = master_edge.first;
//
//      chi_log.Log() << "Master edge identified "
//                    << master_edge.first << "->" << master_edge.second;
//
//      const auto& master_vertex = *mesh.vertices[master_vid];
//
////      vertex_ids.insert(master_edge.first);
//      vertex_ids_set.insert(master_edge.second);
//
//      for (auto& cell_edge_pair : cell_edge_pairs_on_plane)
//        if (cell_edge_pair.first == master_cell_edge_pair.first)
//          continue;
//        else
//        {
//          auto& slave_cell = *cell_edge_pair.first;
//          int   slave_edge_index = cell_edge_pair.second;
//          auto slave_edge = MakeEdgeFromPolygonEdgeIndex(slave_cell,slave_edge_index);
//
//          vertex_ids_set.insert(slave_edge.first);
//          vertex_ids_set.insert(slave_edge.second);
//        }
//
//      std::vector<uint64_t> vertex_ids(vertex_ids_set.begin(),
//                                       vertex_ids_set.end());
//
//      auto VertexComparator = [&mesh,&master_vertex](
//        const uint64_t& v0id, const uint64_t& v1id)
//      {
//        const auto& v0 = *mesh.vertices[v0id];
//        const auto& v1 = *mesh.vertices[v1id];
//
//        double dv0 = (v0 - master_vertex).NormSquare();
//        double dv1 = (v1 - master_vertex).NormSquare();
//
//        return dv0 < dv1;
//      };
//
//      std::sort(vertex_ids.begin(),vertex_ids.end(),VertexComparator);
//
//      chi_log.Log() << "Sorted list of vids relative to " << master_vid <<":";
//      for (auto vid : vertex_ids)
//        chi_log.Log() << vid;
//    }
  }//If more than 2 cells


//  //============================================= Verbose output
//  for (auto cell_ptr : cells_to_add_to_mesh)
//  {
//    chi_log.Log() << "new_cell:";
//    for (auto vid : cell_ptr->vertex_ids)
//      chi_log.Log() << vid;
//  }
//
//  chi_log.Log() << "Done processing cell cut";

}