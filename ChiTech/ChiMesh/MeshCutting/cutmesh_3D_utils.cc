#include "meshcutting.h"

#include <algorithm>
#include <stdexcept>

//###################################################################
/**Checks the quality of a polyhedron against the following.
 *
 * Take each face and form triangles using the face edges and the face
 * centroid. Each triangle must ccw-orientation and have a normal pointing
 * in the same general direction as the average face normal.
 *
 * Now take each triangle and form a tetrahedron with the cell-centroid. The
 * quality requirement is that no inverted tetrahedron face is inverted with
 * respect to the cell-centroid.*/
bool chi_mesh::mesh_cutting::
  CheckPolyhedronQuality(const MeshContinuum &mesh, const chi_mesh::Cell& cell)
{
  const auto& C = cell.centroid;

  for (auto& face : cell.faces)
  {
    const auto& v0 = face.centroid;
    size_t num_edges = face.vertex_ids.size();

    for (int e=0; e<num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids,e);

      const auto& v1 = *mesh.vertices[edge.first];
      const auto& v2 = *mesh.vertices[edge.second];

      auto v01 = v1-v0;
      auto v02 = v2-v0;

      auto n = v01.Cross(v02);
      auto C_tri = (v0+v1+v2)/3.0;

      auto CC = C_tri - C;

      if (CC.Dot(n)<0.0)
        return false;
    }
  }//for face

  return true;
}

//###################################################################
/***/
void chi_mesh::mesh_cutting::
  CutTetrahedron(const std::vector<ECI> &cut_edges,
                 const std::set<uint64_t> &cut_vertices,
                 const Vector3 &plane_point,
                 const Vector3 &plane_normal,
                 MeshContinuum &mesh,
                 chi_mesh::CellPolyhedron &cell)
{
  const std::string fname = __FUNCTION__;
  const auto& p = plane_point;
  const auto& n = plane_normal;

  const size_t num_faces = cell.faces.size();

  /**Utility lambda to check if a vertex is in "cut_vertices" list.*/
  auto VertexIsCut = [&cut_vertices](uint64_t vid)
  {
    auto result = cut_vertices.find(vid);

    if (result != cut_vertices.end())
      return true;

    return false;
  };

  /**Utility function to check if an edge is in the "cut_edges" list.*/
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

  //======================================== Determine number of vertices cut
  size_t num_vertices_cut = 0;

  for (auto vid : cell.vertex_ids)
    if (VertexIsCut(vid))
      ++num_vertices_cut;

  //======================================== Determine number of edges cut
  typedef std::pair<bool,ECI> EdgeCutCon; //Edge Cut Condition
  typedef std::vector<EdgeCutCon> FCC;    //Face Cut Condition per edge

  size_t num_edges_cut = 0;
  std::vector<FCC> CCC(num_faces); //Cell Cut Condition per face
  for (int f=0; f<num_faces; ++f)
  {
    const auto& face = cell.faces[f];

    size_t num_edges = face.vertex_ids.size();
    CCC[f].reserve(num_edges);

    for (int e=0; e<num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids,e);

      auto cut_condition = EdgeIsCut(edge);

      if (cut_condition.first)
      {
        ++num_edges_cut;
        CCC[f].push_back(cut_condition);
      }
      else
        CCC[f].emplace_back(false,ECI());

    }//for e
  }//for faces

  //======================================== Handle case by case
  typedef std::vector<uint64_t>   RawPolygon;
  typedef std::vector<RawPolygon> RawFaces;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 1
  // 3 Edges 0 Verts indicates that a single
  // face is not cut.
  if (num_edges_cut == 3 and num_vertices_cut == 0)
  {
    RawFaces raw_faces;

    //================================= Find base face and cut faces
    int               base_face_index = -1;
    std::vector<int>  cut_faces_indices; cut_faces_indices.reserve(3);
    for (int f=0; f<num_faces; ++f)
    {
      bool no_edges_cut = true;
      size_t num_edges = cell.faces[f].vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        if (CCC[f][e].first) {no_edges_cut = false; break;}

      if (no_edges_cut) {base_face_index = f; break;}
      else              {cut_faces_indices.push_back(f);}
    }//for f, find base face
    if (base_face_index < 0)
      throw std::logic_error(fname + ": Case 1 base_face_index < 0.");

    raw_faces.push_back(cell.faces[base_face_index].vertex_ids);

    //================================= Process cut faces
    for (auto cf_index : cut_faces_indices)
    {
      auto& face = cell.faces[cf_index];
      size_t num_edges = face.vertex_ids.size();

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<num_edges; ++e)
      {
        int ep1 = (e<(num_edges-1))? e+1 : 0;

        bool curr_edge_cut = CCC[cf_index][e].first;
        bool next_edge_cut = CCC[cf_index][ep1].first;

        ECI curr_eci = CCC[cf_index][e].second;

        if (not curr_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
        else if (curr_edge_cut and next_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.cut_point_id);
        }
        else if (curr_edge_cut and (not next_edge_cut))
        {
          raw_face.push_back(curr_eci.cut_point_id);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
      }//for e
      raw_faces.push_back(raw_face);
    }//for cf_index


  }//Case 1
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 2
  else if (num_edges_cut == 2 and num_vertices_cut == 1)
  {

  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 3
  else if (num_edges_cut == 1 and num_vertices_cut == 2)
  {

  }
  else
    throw std::logic_error(fname + ": Unsupported cut case.");

}