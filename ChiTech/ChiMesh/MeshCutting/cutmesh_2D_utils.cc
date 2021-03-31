#include "meshcutting.h"

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
/***/
void chi_mesh::mesh_cutting::
  SplitConcavePolygons(MeshContinuum &mesh,
                       std::vector<chi_mesh::Cell*> &cell_list)
{
  const chi_mesh::Vector3 khat(0.0,0.0,1.0);

  //======================================== Lambda to find a vertex
  //                                         in and edge-list
  auto EdgeIndexWithVertexAtEnd = [](
    const std::vector<Edge>& edge_list,
    const uint64_t vertex_id)
  {
    int e=0;
    for (auto& edge : edge_list)
    {
      if (edge.second == vertex_id)
        return e;
      ++e;
    }
    return -1;
  };

  //======================================== Make copy of cell-list
  std::vector<chi_mesh::Cell*> original_cell_list = cell_list;

  //======================================== Loop over cells
  for (auto& cell_ptr : original_cell_list)
  {
    auto&  cell      = *(chi_mesh::CellPolygon*)cell_ptr;
    size_t num_edges = cell.vertex_ids.size();
    size_t num_verts = num_edges;

    //================================= Make edges
    std::vector<Edge> cell_edges(num_edges);
    for (int e=0; e<num_edges; ++e)
      cell_edges[e] = MakeEdgeFromPolygonEdgeIndex(cell,e);

    //================================= Loop over vertices and flag
    //                                  concave inflections
    std::vector<bool> vertex_at_concavity(num_verts,false);
    for (int v=0; v<num_verts; v++)
    {
      int e0 = EdgeIndexWithVertexAtEnd(cell_edges,cell.vertex_ids[v]);
      int e1 = (e0 < (num_edges-1))? e0+1 : 0;

      const auto& edge0 = cell_edges[e0];
      const auto& edge1 = cell_edges[e1];

      const auto& v0 = *mesh.vertices[edge0.first];
      const auto& v1 = *mesh.vertices[edge0.second];
      const auto& v2 = *mesh.vertices[edge1.second];

      auto v01 = v1-v0;
      auto v12 = v2-v1;

      if (v01.Cross(v12).Dot(khat)<0.0)
        vertex_at_concavity[v] = true;
    }//for v


  }//for cell
}