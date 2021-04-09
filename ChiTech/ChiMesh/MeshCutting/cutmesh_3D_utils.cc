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
/**Defines a polyhedron based only on its faces.*/
void chi_mesh::mesh_cutting::
  PopulatePolyhedronFromFaces(const MeshContinuum &mesh,
                              const std::vector<std::vector<uint64_t>> &raw_faces,
                              chi_mesh::Cell &cell)
{
  const size_t num_raw_faces = raw_faces.size();
  std::vector<uint64_t> cell_vertex_ids;

  /**Lamda to add a vertex to a list.*/
  auto AddIDToCellVertexIDs = [&cell_vertex_ids](uint64_t new_id)
  {
    auto find_result = std::find(cell_vertex_ids.begin(),
                                 cell_vertex_ids.end(),
                                 new_id);

    if (find_result == cell_vertex_ids.end())
      cell_vertex_ids.push_back(new_id);
  };

  //======================================== Build faces
  cell.faces.clear();
  for (auto& raw_face : raw_faces)
  {
    chi_mesh::Vector3 face_centroid;
    for (uint64_t vid : raw_face)
    {
      face_centroid += *mesh.vertices[vid];
      AddIDToCellVertexIDs(vid);
    }
    face_centroid /= (double)raw_face.size();

    chi_mesh::CellFace new_face;
    new_face.centroid = face_centroid;
    new_face.vertex_ids = raw_face;

    cell.faces.push_back(std::move(new_face));
  }//for raw face

  //======================================== Compute cell centroid
  chi_mesh::Vector3 cell_centroid;
  for (uint64_t vid : cell_vertex_ids)
    cell_centroid += *mesh.vertices[vid];
  cell_centroid /= (double)cell_vertex_ids.size();

  cell.centroid = cell_centroid;
  cell.vertex_ids = cell_vertex_ids;
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

  //======================================== Determine number of edges cut
  //                                         and number of vertices cut
  size_t num_edges_cut = cut_edges.size();
  size_t num_vertices_cut = cut_vertices.size();

  //======================================== Determine number of edges cut
  typedef std::pair<bool,ECI> EdgeCutCon; //Edge Cut Condition
  typedef std::vector<EdgeCutCon> FCC;    //Face Cut Condition per edge

  std::vector<FCC> CCC(num_faces); //Cell Cut Condition per face
  for (int f=0; f<num_faces; ++f)
  {
    const auto& face = cell.faces[f];

    size_t num_edges = face.vertex_ids.size();
    CCC[f].resize(num_edges,EdgeCutCon(false,ECI()));

    for (int e=0; e<num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids,e);
      auto cut_condition = EdgeIsCut(edge);

      if (cut_condition.first) CCC[f][e] = cut_condition;
      else                     CCC[f][e] = std::make_pair(false,ECI());
    }//for e
  }//for faces

  //======================================== Handle case by case
  typedef std::vector<uint64_t>   RawPolygon;
  typedef std::vector<RawPolygon> RawFaces;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 1
  // 3 Edges 0 Verts indicates that a single
  // face is not cut. Since we are cutting a
  // Tet this can only mean that we created
  // a polyhedron and a tet sliver with a
  // triangle as an interface-face.
  // Nomenclature:
  //   IF = interface-face
  //   BF = base-face
  //   CF = cut-face
  if (num_edges_cut == 3 and num_vertices_cut == 0)
  {
    RawFaces raw_faces_polyhedron;
    RawFaces raw_faces_slivertet;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);

    if (IF_vids.size() != 3)
      throw std::logic_error(fname + ": Case 1 IF_vids.size() != 3.");

    //================================= Find base face and cut faces
    // The base face is the face that
    // was not cut. All other faces
    // are cut faces
    int               BF_id = -1;
    std::vector<int>  CF_ids;
    CF_ids.reserve(3);
    for (int f=0; f<num_faces; ++f)
    {
      bool no_edges_cut = true;
      size_t num_edges = cell.faces[f].vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        if (CCC[f][e].first) {no_edges_cut = false; break;}

      if (no_edges_cut) { BF_id = f;}
      else              {CF_ids.push_back(f);}
    }//for f, find base face
    if (BF_id < 0)
      throw std::logic_error(fname + ": Case 1 base_face_index < 0.");

    const auto& BF = cell.faces[BF_id]; //Base-face

    //================================= Establish order of vertices
    //                                  for polyhedron and sliver-Tet

    // Compute interface-face centroid (IFC)
    chi_mesh::Vector3 IFC;
    for (auto ivid : IF_vids)
      IFC += *mesh.vertices[ivid];
    IFC /= (double)IF_vids.size();

    // Form vector from base-face centroid to
    // interface-face centroid
    auto BFC_IFC = IFC - BF.centroid;

    // Take existing interface-face vertex ordering
    // and build a normal for the interface-face
    const auto& IFv0 = *mesh.vertices[IF_vids[0]];
    const auto& IFv1 = *mesh.vertices[IF_vids[1]];
    const auto& IFv2 = *mesh.vertices[IF_vids[2]];

    auto IFv01 = IFv1 - IFv0;
    auto IFv02 = IFv2 - IFv0;

    auto IFn = IFv01.Cross(IFv02);

    // Now check the existing orientation
    // with respect to the base-face. If
    // the orientation is wrong, flip it.
    if (IFn.Dot(BFC_IFC) < 0.0)
      IF_vids = {IF_vids[0],IF_vids[2],IF_vids[1]};

    std::vector<uint64_t> IF_vids_reverse(IF_vids.rbegin(),IF_vids.rend());

    //================================= Process cut faces to form faces
    //                                  for the polyhedron
    // This process also identifies the
    // the dangling vertex that form
    // part of the sliver Tet
    uint64_t dangling_vertex_id = 0;
    for (int CF_id : CF_ids)
    {
      auto& face = cell.faces[CF_id];
      size_t num_edges = face.vertex_ids.size();

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<num_edges; ++e)
      {
        int ep1 = (e<(num_edges-1))? e+1 : 0;

        bool curr_edge_cut = CCC[CF_id][e]  .first;
        bool next_edge_cut = CCC[CF_id][ep1].first;

        ECI curr_eci = CCC[CF_id][e].second;

        if (not curr_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
        else if (curr_edge_cut and next_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.cut_point_id);
          dangling_vertex_id = curr_eci.vertex_ids.second;
        }
        else if (curr_edge_cut and (not next_edge_cut))
        {
          raw_face.push_back(curr_eci.cut_point_id);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
      }//for e
      raw_faces_polyhedron.push_back(raw_face);
    }//for cf_index

    //================================= Add base-face and interface-face
    //                                  to polyhedron raw faces
    raw_faces_polyhedron.push_back(cell.faces[BF_id].vertex_ids);
    raw_faces_polyhedron.push_back(IF_vids);

    //================================= Build the faces of the sliver-tet
    raw_faces_slivertet.push_back(IF_vids_reverse);

    size_t IF_num_edges = IF_vids_reverse.size();
    for (int e=0; e<IF_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(IF_vids_reverse,e);

      RawPolygon raw_face = {edge.first,edge.second,dangling_vertex_id};
      raw_faces_slivertet.emplace_back(raw_face);
    }//for e

    //================================= Complete cell definitions
    PopulatePolyhedronFromFaces(mesh,raw_faces_polyhedron,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_slivertet,*new_cell);

    mesh.cells.push_back(new_cell);
  }//end Case 1
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 2
  // 2 Edges 1 Verts indicates again 1 uncut
  // face. The cut will result in a polyhedron
  // and tet sliver.
  else if (num_edges_cut == 2 and num_vertices_cut == 1)
  {
    RawFaces raw_faces_polyhedron;
    RawFaces raw_faces_slivertet;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);
    for (auto& cut_vertex : cut_vertices)
      IF_vids.push_back(cut_vertex);

    if (IF_vids.size() != 3)
      throw std::logic_error(fname + ": Case 2 IF_vids.size() != 3.");

    //================================= Find base face and cut faces
    // The base face is the face that
    // was not cut. All other faces
    // are cut faces
    int               BF_id = -1;
    std::vector<int>  CF_ids;
    CF_ids.reserve(3);
    for (int f=0; f<num_faces; ++f)
    {
      bool no_edges_cut = true;
      size_t num_edges = cell.faces[f].vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        if (CCC[f][e].first) {no_edges_cut = false; break;}

      if (no_edges_cut) { BF_id = f;}
      else              {CF_ids.push_back(f);}
    }//for f, find base face
    if (BF_id < 0)
      throw std::logic_error(fname + ": Case 2 base_face_index < 0.");

    const auto& BF = cell.faces[BF_id]; //Base-face

    //================================= Establish order of vertices
    //                                  for polyhedron and sliver-Tet

    // Compute interface-face centroid (IFC)
    chi_mesh::Vector3 IFC;
    for (auto ivid : IF_vids)
      IFC += *mesh.vertices[ivid];
    IFC /= (double)IF_vids.size();

    // Form vector from base-face centroid to
    // interface-face centroid
    auto BFC_IFC = IFC - BF.centroid;

    // Take existing interface-face vertex ordering
    // and build a normal for the interface-face
    const auto& IFv0 = *mesh.vertices[IF_vids[0]];
    const auto& IFv1 = *mesh.vertices[IF_vids[1]];
    const auto& IFv2 = *mesh.vertices[IF_vids[2]];

    auto IFv01 = IFv1 - IFv0;
    auto IFv02 = IFv2 - IFv0;

    auto IFn = IFv01.Cross(IFv02);

    // Now check the existing orientation
    // with respect to the base-face. If
    // the orientation is wrong, flip it.
    if (IFn.Dot(BFC_IFC) < 0.0)
      IF_vids = {IF_vids[0],IF_vids[2],IF_vids[1]};

    std::vector<uint64_t> IF_vids_reverse(IF_vids.rbegin(),IF_vids.rend());

    //================================= Process cut faces to form faces
    //                                  for the polyhedron
    // This process also identifies the
    // the dangling vertex that form
    // part of the sliver Tet.
    // - 2 cut-faces. 1 edge is cut, 1 vertex is cut
    // - 1 cut-face. 2 edges are cut.
    uint64_t dangling_vertex_id = 0;
    for (int CF_id : CF_ids)
    {
      auto& face = cell.faces[CF_id];
      size_t num_edges = face.vertex_ids.size();

      size_t face_num_cut_edges = 0;
      for (int e=0; e<num_edges; ++e)
        if (CCC[CF_id][e].first)
          ++face_num_cut_edges;

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<num_edges; ++e)
      {
        int ep1 = (e<(num_edges-1))? e+1 : 0;

        bool curr_edge_cut = CCC[CF_id][e].first;

        ECI curr_eci = CCC[CF_id][e].second;

        if (face_num_cut_edges == 1)
        {
          bool end_vertex_cut = VertexIsCut(ep1);

          if ((not curr_edge_cut) and (not end_vertex_cut))
          {
            raw_face.push_back(curr_eci.vertex_ids.first);
            raw_face.push_back(curr_eci.vertex_ids.second);
          }
          else if (curr_edge_cut)
          {
            raw_face.push_back(curr_eci.vertex_ids.first);
            raw_face.push_back(curr_eci.cut_point_id);
            dangling_vertex_id = curr_eci.vertex_ids.second;
          }
          // The third edge is discarded
        }
        else if (face_num_cut_edges == 2)
        {
          bool next_edge_cut = CCC[CF_id][ep1].first;

          if (not curr_edge_cut)
          {
            raw_face.push_back(curr_eci.vertex_ids.first);
            raw_face.push_back(curr_eci.vertex_ids.second);
          }
          else if (curr_edge_cut and next_edge_cut)
          {
            raw_face.push_back(curr_eci.vertex_ids.first);
            raw_face.push_back(curr_eci.cut_point_id);
            dangling_vertex_id = curr_eci.vertex_ids.second;
          }
          else if (curr_edge_cut and (not next_edge_cut))
          {
            raw_face.push_back(curr_eci.cut_point_id);
            raw_face.push_back(curr_eci.vertex_ids.second);
          }
        }
        else
          throw std::logic_error(fname + ": Case 2 face_num_cut_edges != 1or2.");

      }//for e
      raw_faces_polyhedron.push_back(raw_face);
    }//for cf_index

    //================================= Add base-face and interface-face
    //                                  to polyhedron raw faces
    raw_faces_polyhedron.push_back(cell.faces[BF_id].vertex_ids);
    raw_faces_polyhedron.push_back(IF_vids);

    //================================= Build the faces of the sliver-tet
    raw_faces_slivertet.push_back(IF_vids_reverse);

    size_t IF_num_edges = IF_vids_reverse.size();
    for (int e=0; e<IF_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(IF_vids_reverse,e);

      RawPolygon raw_face = {edge.first,edge.second,dangling_vertex_id};
      raw_faces_slivertet.emplace_back(raw_face);
    }//for e

    //================================= Complete cell definitions
    PopulatePolyhedronFromFaces(mesh,raw_faces_polyhedron,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_slivertet,*new_cell);

    mesh.cells.push_back(new_cell);
  }//end Case 2
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 3
  // This type of cut results in two tets
  else if (num_edges_cut == 1 and num_vertices_cut == 2)
  {
    RawFaces raw_faces_tetA;
    RawFaces raw_faces_tetB;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);
    for (auto& cut_vertex : cut_vertices)
      IF_vids.push_back(cut_vertex);

    if (IF_vids.size() != 3)
      throw std::logic_error(fname + ": Case 3 IF_vids.size() != 3.");

    //================================= Find base face and cut faces
    // The base face is the face that
    // was not cut. All other faces
    // are cut faces
    int               BF_id = -1;
    std::vector<int>  CF_ids;
    CF_ids.reserve(3);
    for (int f=0; f<num_faces; ++f)
    {
      bool no_edges_cut = true;
      size_t num_edges = cell.faces[f].vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
        if (CCC[f][e].first) {no_edges_cut = false; break;}

      if (no_edges_cut) { BF_id = f;}
      else              {CF_ids.push_back(f);}
    }//for f, find base face
    if (BF_id < 0)
      throw std::logic_error(fname + ": Case 3 base_face_index < 0.");

    const auto& BF = cell.faces[BF_id]; //Base-face

    //================================= Establish order of vertices
    //                                  for polyhedron and sliver-Tet

    // Compute interface-face centroid (IFC)
    chi_mesh::Vector3 IFC;
    for (auto ivid : IF_vids)
      IFC += *mesh.vertices[ivid];
    IFC /= (double)IF_vids.size();

    // Form vector from base-face centroid to
    // interface-face centroid
    auto BFC_IFC = IFC - BF.centroid;

    // Take existing interface-face vertex ordering
    // and build a normal for the interface-face
    const auto& IFv0 = *mesh.vertices[IF_vids[0]];
    const auto& IFv1 = *mesh.vertices[IF_vids[1]];
    const auto& IFv2 = *mesh.vertices[IF_vids[2]];

    auto IFv01 = IFv1 - IFv0;
    auto IFv02 = IFv2 - IFv0;

    auto IFn = IFv01.Cross(IFv02);

    // Now check the existing orientation
    // with respect to the base-face. If
    // the orientation is wrong, flip it.
    if (IFn.Dot(BFC_IFC) < 0.0)
      IF_vids = {IF_vids[0],IF_vids[2],IF_vids[1]};

    std::vector<uint64_t> IF_vids_reverse(IF_vids.rbegin(),IF_vids.rend());

    //================================= Process cut faces to form faces
    //                                  for the polyhedron
    // This process also identifies the
    // the dangling vertex that form
    // part of the sliver Tet
    uint64_t dangling_vertex_id = 0;
    for (int CF_id : CF_ids)
    {
      auto& face = cell.faces[CF_id];
      size_t num_edges = face.vertex_ids.size();

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<num_edges; ++e)
      {
        int ep1 = (e<(num_edges-1))? e+1 : 0;

        bool curr_edge_cut = CCC[CF_id][e]  .first;
        bool next_edge_cut = CCC[CF_id][ep1].first;

        ECI curr_eci = CCC[CF_id][e].second;

        if (not curr_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
        else if (curr_edge_cut and next_edge_cut)
        {
          raw_face.push_back(curr_eci.vertex_ids.first);
          raw_face.push_back(curr_eci.cut_point_id);
          dangling_vertex_id = curr_eci.vertex_ids.second;
        }
        else if (curr_edge_cut and (not next_edge_cut))
        {
          raw_face.push_back(curr_eci.cut_point_id);
          raw_face.push_back(curr_eci.vertex_ids.second);
        }
      }//for e
      raw_faces_tetA.push_back(raw_face);
    }//for cf_index

    //================================= Add base-face and interface-face
    //                                  to polyhedron raw faces
    raw_faces_tetA.push_back(cell.faces[BF_id].vertex_ids);
    raw_faces_tetA.push_back(IF_vids);

    //================================= Build the faces of the sliver-tet
    raw_faces_tetB.push_back(IF_vids_reverse);

    size_t IF_num_edges = IF_vids_reverse.size();
    for (int e=0; e<IF_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(IF_vids_reverse,e);

      RawPolygon raw_face = {edge.first,edge.second,dangling_vertex_id};
      raw_faces_tetB.emplace_back(raw_face);
    }//for e

    //================================= Complete cell definitions
    PopulatePolyhedronFromFaces(mesh,raw_faces_tetA,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_tetB,*new_cell);

    mesh.cells.push_back(new_cell);
  }
  else
    throw std::logic_error(fname + ": Unsupported cut case.");

}