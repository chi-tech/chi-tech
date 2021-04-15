#include "meshcutting.h"

#include <algorithm>
#include <stdexcept>

#include "chi_log.h"

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
  CutTetrahedron(const std::vector<ECI> &global_cut_edges,
                 const std::set<uint64_t> &global_cut_vertices,
                 const Vector3 &plane_point,
                 const Vector3 &plane_normal,
                 MeshContinuum &mesh,
                 chi_mesh::CellPolyhedron &cell)
{
  const std::string fname = __FUNCTION__;
  const auto& p = plane_point;
  const auto& n = plane_normal;

  const size_t num_faces = cell.faces.size();

  ChiLog& chi_log = ChiLog::GetInstance();

  /**Utility lambda to check if a vertex is in "cut_vertices" list.*/
  auto VertexIsCut = [&global_cut_vertices](uint64_t vid)
  {
    auto result = global_cut_vertices.find(vid);

    if (result != global_cut_vertices.end())
      return true;

    return false;
  };

  /**Utility function to check if an edge is in the "cut_edges" list.*/
  auto EdgeIsCut = [&global_cut_edges](const Edge& edge)
  {
    Edge edge_set(std::min(edge.first,edge.second),
                  std::max(edge.first,edge.second));

    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto result = std::find_if(global_cut_edges.begin(),
                               global_cut_edges.end(),
                               std::bind(Comparator,Arg1,edge_set));

    if (result != global_cut_edges.end())
      return std::make_pair(true,*result);

    return std::make_pair(false,*result);
  };

  auto EdgesHaveSameVerts = [](const Edge& edgeA, const Edge& edgeB)
  {
    bool criss_cross = ((edgeA.first == edgeB.second) and
                        (edgeA.second == edgeB.first));
    bool straight = (edgeA == edgeB);

    return (straight or criss_cross);
  };

  bool verbose = true;
  //======================================== Verbose output
  if (verbose)
  {
    chi_log.Log() << "Cell:";
    for (uint64_t vid : cell.vertex_ids)
      chi_log.Log() << vid << " " << mesh.vertices[vid]->PrintS();
  }

  //======================================== Build cell unique edges
  std::set<Edge> edges_set;
  for (const auto& face : cell.faces)
  {
    const size_t face_num_edges = face.vertex_ids.size();

    for (int e=0; e<face_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids,e);
      edges_set.insert(std::make_pair(std::min(edge.first,edge.second),
                                      std::max(edge.first,edge.second)));
    }
  }//for face

  //======================================== Determine number of edges cut
  //                                         and number of vertices cut
  std::vector<ECI> cell_cut_edges;
  for (const auto& edge : edges_set)
  {
    auto edge_cut_condition = EdgeIsCut(edge);
    if (edge_cut_condition.first)
      cell_cut_edges.push_back(edge_cut_condition.second);
  }
  size_t num_edges_cut = cell_cut_edges.size();

  std::vector<uint64_t> cell_cut_vertices;
  for (uint64_t vid : cell.vertex_ids)
    if (VertexIsCut(vid)) cell_cut_vertices.push_back(vid);
  size_t num_vertices_cut = cell_cut_vertices.size();


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
      auto  edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids,e);
      auto  cut_condition = EdgeIsCut(edge);
      auto& eci = cut_condition.second;

      if (cut_condition.first)
        CCC[f][e] = std::make_pair(true,ECI(edge,eci.cut_point_id));
      else
        CCC[f][e] = std::make_pair(false,ECI(edge,0));
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
    if (verbose)
      chi_log.Log(LOG_0VERBOSE_1) << "Cell split case 1";
    RawFaces raw_faces_polyhedron;
    RawFaces raw_faces_slivertet;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cell_cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);
    for (auto& cut_vertex : cell_cut_vertices)
      IF_vids.push_back(cut_vertex);

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

    //================================= Establish order of interface-face
    //                                  vertices wrt polyhedron and sliver-Tet
    // This is fairly easy because we
    // know the interface-face is a
    // triangle.

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

    if (verbose)
    {
      chi_log.Log() << "IF:";
      for (uint64_t vid : IF_vids)
        chi_log.Log() << vid << " " << mesh.vertices[vid]->PrintS();
      chi_log.Log() << "IFrev:";
      for (uint64_t vid : IF_vids_reverse)
        chi_log.Log() << vid << " " << mesh.vertices[vid]->PrintS();
    }

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

        if (curr_edge_cut and next_edge_cut)
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

      RawPolygon raw_face = {edge.second,edge.first,dangling_vertex_id};
      raw_faces_slivertet.emplace_back(raw_face);
    }//for e

    //================================= Complete cell definitions
    PopulatePolyhedronFromFaces(mesh,raw_faces_polyhedron,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_slivertet,*new_cell);

    mesh.cells.push_back(new_cell);

    bool checkA = CheckPolyhedronQuality(mesh, cell);
    bool checkB = CheckPolyhedronQuality(mesh, *new_cell);

    if (verbose)
    {
      chi_log.Log() << "Cell A.";
      for (auto& face : cell.faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
      chi_log.Log() << "Cell B.";
      for (auto& face : new_cell->faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
    }

    if (not checkA)
      throw std::logic_error(fname + ": Cell quality check A failed");
    if (not checkB)
      throw std::logic_error(fname + ": Cell quality check B failed");
  }//end Case 1
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 2
  // 2 Edges 1 Verts indicates again 1 uncut
  // face. The cut will result in a polyhedron
  // and tet sliver.
  else if (num_edges_cut == 2 and num_vertices_cut == 1)
  {
    chi_log.Log(LOG_0VERBOSE_1) << "Cell split case 2";
    RawFaces raw_faces_polyhedron;
    RawFaces raw_faces_slivertet;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cell_cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);
    for (auto& cut_vertex : cell_cut_vertices)
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

    if (verbose)
    {
      chi_log.Log() << "IF:";
      for (uint64_t vid : IF_vids)
        chi_log.Log() << vid << " " << mesh.vertices[vid]->PrintS();
      chi_log.Log() << "IFrev:";
      for (uint64_t vid : IF_vids_reverse)
        chi_log.Log() << vid << " " << mesh.vertices[vid]->PrintS();
    }

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
          bool end_vertex_cut = VertexIsCut(curr_eci.vertex_ids.second);

          if ((not curr_edge_cut) and (not end_vertex_cut))
          {
            raw_face.push_back(curr_eci.vertex_ids.first);
//            raw_face.push_back(curr_eci.vertex_ids.second);
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
//            raw_face.push_back(curr_eci.vertex_ids.first);
            raw_face.push_back(curr_eci.vertex_ids.second);
          }
          else if (curr_edge_cut and next_edge_cut)
          {
//            raw_face.push_back(curr_eci.vertex_ids.first);
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

    bool checkA = CheckPolyhedronQuality(mesh, cell);
    bool checkB = CheckPolyhedronQuality(mesh, *new_cell);

    if (verbose)
    {
      chi_log.Log() << "Cell A.";
      for (auto& face : cell.faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
      chi_log.Log() << "Cell B.";
      for (auto& face : new_cell->faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
    }

    if (not checkA)
      throw std::logic_error(fname + ": Cell quality check A failed");
    if (not checkB)
      throw std::logic_error(fname + ": Cell quality check B failed");
  }//end Case 2
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 3
  // This type of cut results in two tets
  else if (num_edges_cut == 1 and num_vertices_cut == 2)
  {
    chi_log.Log(LOG_0VERBOSE_1) << "Cell split case 3";
    RawFaces raw_faces_tetA;
    RawFaces raw_faces_tetB;

    //================================= Build interface face
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(3);
    for (auto& cut_edge : cell_cut_edges)
      IF_vids.push_back(cut_edge.cut_point_id);
    for (auto& cut_vertex : cell_cut_vertices)
      IF_vids.push_back(cut_vertex);

    if (IF_vids.size() != 3)
      throw std::logic_error(fname + ": Case 3 IF_vids.size() != 3.");

    //================================= Find dangling vertices
    std::vector<uint64_t> dangling_vertex_ids;
    dangling_vertex_ids.reserve(2);
    for (uint64_t vid : cell.vertex_ids)
      if (!VertexIsCut(vid)) dangling_vertex_ids.push_back(vid);

    if (dangling_vertex_ids.size() != 2)
      throw std::logic_error(fname + ": Case 3 dangling_vertices.size() != 2.");

    //================================= Establish order of vertices
    //                                  for polyhedron and sliver-Tet

    // Compute interface-face centroid (IFC)
    chi_mesh::Vector3 IFC;
    for (auto ivid : IF_vids)
      IFC += *mesh.vertices[ivid];
    IFC /= (double)IF_vids.size();

    // Form vector from DV0 to
    // interface-face centroid
    const auto& DV0 = *mesh.vertices[dangling_vertex_ids[0]];
    auto BFC_IFC = IFC - DV0;

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

    //================================= Add interface-faces to tets
    raw_faces_tetA.push_back(IF_vids);
    raw_faces_tetB.push_back(IF_vids_reverse);

    //================================= Build the other faces of tet-A
    size_t IF_num_edges = IF_vids.size();
    for (int e=0; e<IF_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(IF_vids,e);

      RawPolygon raw_face = {edge.second,edge.first,dangling_vertex_ids[0]};
      raw_faces_tetA.emplace_back(raw_face);
    }//for e

    //================================= Build the other faces of tet-B;
    for (int e=0; e<IF_num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(IF_vids_reverse,e);

      RawPolygon raw_face = {edge.second,edge.first,dangling_vertex_ids[1]};
      raw_faces_tetB.emplace_back(raw_face);
    }//for e

    //================================= Complete cell definitions
    PopulatePolyhedronFromFaces(mesh,raw_faces_tetA,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_tetB,*new_cell);

    mesh.cells.push_back(new_cell);

    bool checkA = CheckPolyhedronQuality(mesh, cell);
    bool checkB = CheckPolyhedronQuality(mesh, *new_cell);

    if (not checkA)
      throw std::logic_error(fname + ": Cell quality check A failed");
    if (not checkB)
      throw std::logic_error(fname + ": Cell quality check B failed");
  }//end Case 3
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 4
  // This results in 2 polyhedra.
  // All faces are cut at 2 points.
  else if (num_edges_cut == 4 and num_vertices_cut == 0)
  {
    chi_log.Log(LOG_0VERBOSE_1) << "Cell split case 4";

    RawFaces raw_faces_polyhedronA;
    RawFaces raw_faces_polyhedronB;

    std::array<uint64_t,2> dangling_verticesA = {0,0};
    std::array<uint64_t,2> dangling_verticesB = {0,0};

    /**Lamda to check if a vertex is a dangling vertex.*/
    auto IsVertexDanglingVertex = [](const std::array<uint64_t,2> dang_verts,
                                     uint64_t vid)
    {
      if (vid == dang_verts[0] or vid == dang_verts[1])
        return true;

      return false;
    };

    //================================= Identify dangling vertices for
    //                                  polyhedra A and B
    size_t num_dv_A = 0;
    size_t num_dv_B = 0;
    for (uint64_t vid : cell.vertex_ids)
    {
      const auto& v = *mesh.vertices[vid];

      if ((v-p).Dot(n) < 0.0) dangling_verticesA[num_dv_A++] = vid;
      else                    dangling_verticesB[num_dv_B++] = vid;
    }

    if (num_dv_A != 2 or num_dv_B != 2)
      throw std::logic_error(fname + ": Case 4 num_dv_A != 2 or num_dv_B != 2.");

    //================================= Compute dangling verticesA centroids
    chi_mesh::Vector3 DVC_A = 0.5*(*mesh.vertices[dangling_verticesA[0]] +
                                   *mesh.vertices[dangling_verticesA[1]]);

    //================================= Build interface-face edges (unordered)
    std::vector<Edge> IF_edges_unordered;
    IF_edges_unordered.reserve(4);
    for (int f=0; f<num_faces; ++f)
    {
      std::array<uint64_t,2> cut_point_ids = {0,0};
      size_t num_cut_points=0;
      size_t num_edges = cell.faces[f].vertex_ids.size();
      for (int e=0; e<num_edges; ++e)
      {
        bool edge_is_cut = CCC[f][e].first;
        ECI& edge_cut_info = CCC[f][e].second;

        if (edge_is_cut)
          cut_point_ids[num_cut_points++] = edge_cut_info.cut_point_id;

        if (num_cut_points>1) break;
      }//for e

      if (num_cut_points != 2)
        throw std::logic_error(fname + ": Case 4 num_cut_points != 2");

      IF_edges_unordered.emplace_back(cut_point_ids[0], cut_point_ids[1]);
    }//for f

    //================================= Order the edges
    std::vector<Edge> IF_edges_ordered;
    IF_edges_ordered.reserve(4);
    IF_edges_ordered.push_back(IF_edges_unordered.back());
    IF_edges_unordered.pop_back();

    for (int i=0; i<3; ++i)
    {
      auto& cur_edge = IF_edges_ordered.back();
      for (auto& other_edge : IF_edges_unordered)
        if (not EdgesHaveSameVerts(cur_edge,other_edge))
        {
          if (other_edge.first == cur_edge.second)
          {IF_edges_ordered.push_back(other_edge); break;}
          if (other_edge.second == cur_edge.second)
          {IF_edges_ordered.emplace_back(other_edge.second,other_edge.first);break;}
        }
    }//for i

    //================================= Determine and correct the interface-face
    //                                  orientation wrt to dangling-vertices A
    // Build vertices
    std::vector<uint64_t> IF_vids;
    IF_vids.reserve(4);
    for (auto& edge : IF_edges_ordered)
      IF_vids.push_back(edge.first);

    // Compute interface-face centroid (IFC)
    chi_mesh::Vector3 IFC;
    for (auto ivid : IF_vids)
      IFC += *mesh.vertices[ivid];
    IFC /= (double)IF_vids.size();

    // Form a vector from dangling vertex centroid
    // to interface-face centroid
    auto DVC_IFC = IFC - DVC_A;

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
    if (IFn.Dot(DVC_IFC) < 0.0)
      IF_vids = {IF_vids[0],IF_vids[3],IF_vids[2],IF_vids[1]};

    std::vector<uint64_t> IF_vids_reverse(IF_vids.rbegin(),IF_vids.rend());

    //================================= Process faces and extract face for
    //                                  polyhedron A only
    raw_faces_polyhedronA.push_back(IF_vids);
    for (int f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      size_t face_num_edges = face.vertex_ids.size();

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<face_num_edges; ++e)
      {
        bool edge_is_cut = CCC[f][e].first;
        auto& eci = CCC[f][e].second;
        auto& edge = eci.vertex_ids;

        uint64_t v0_id = edge.first;
        uint64_t v1_id = edge.second;

        bool v0_dangling = IsVertexDanglingVertex(dangling_verticesA,v0_id);
        bool v1_dangling = IsVertexDanglingVertex(dangling_verticesA,v1_id);

        if (not edge_is_cut)
        {
          if (v0_dangling and v1_dangling)
            raw_face.push_back(v0_id);
        }
        else
        {
          if (v0_dangling and (not v1_dangling))
          {
            raw_face.push_back(v0_id);
            raw_face.push_back(eci.cut_point_id);
          }

          if ((not v0_dangling) and v1_dangling)
            raw_face.push_back(eci.cut_point_id);
        }

      }//for e
      raw_faces_polyhedronA.push_back(raw_face);
    }//for f

    //================================= Process faces and extract face for
    //                                  polyhedron B only
    raw_faces_polyhedronB.push_back(IF_vids_reverse);
    for (int f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      size_t face_num_edges = face.vertex_ids.size();

      RawPolygon raw_face;
      raw_face.reserve(4);
      for (int e=0; e<face_num_edges; ++e)
      {
        bool edge_is_cut = CCC[f][e].first;
        auto& eci = CCC[f][e].second;
        auto& edge = eci.vertex_ids;

        uint64_t v0_id = edge.first;
        uint64_t v1_id = edge.second;

        bool v0_dangling = IsVertexDanglingVertex(dangling_verticesB,v0_id);
        bool v1_dangling = IsVertexDanglingVertex(dangling_verticesB,v1_id);

        if (not edge_is_cut)
        {
          if (v0_dangling and v1_dangling)
            raw_face.push_back(v0_id);
        }
        else
        {
          if (v0_dangling and (not v1_dangling))
          {
            raw_face.push_back(v0_id);
            raw_face.push_back(eci.cut_point_id);
          }

          if ((not v0_dangling) and v1_dangling)
            raw_face.push_back(eci.cut_point_id);
        }

      }//for e
      raw_faces_polyhedronB.push_back(raw_face);
    }//for f

    PopulatePolyhedronFromFaces(mesh,raw_faces_polyhedronA,cell);

    auto new_cell = new chi_mesh::CellPolyhedron;
    PopulatePolyhedronFromFaces(mesh,raw_faces_polyhedronB,*new_cell);

    mesh.cells.push_back(new_cell);

    bool checkA = CheckPolyhedronQuality(mesh, cell);
    bool checkB = CheckPolyhedronQuality(mesh, *new_cell);

    if (verbose)
    {
      chi_log.Log() << "Cell A.";
      for (auto& face : cell.faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
      chi_log.Log() << "Cell B.";
      for (auto& face : new_cell->faces)
      {
        chi_log.Log() << "Face:";
        for (auto vid : face.vertex_ids)
          chi_log.Log() << vid;
      }
    }

    if (not checkA)
      throw std::logic_error(fname + ": Cell quality check A failed");
    if (not checkB)
      throw std::logic_error(fname + ": Cell quality check B failed");
  }//end Case 4
  else
    throw std::logic_error(fname + ": Unsupported cut case. "
                           "num_edges_cut=" + std::to_string(num_edges_cut) +
                           " num_vertices_cut=" + std::to_string(num_vertices_cut));

}