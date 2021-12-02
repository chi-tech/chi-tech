#include "volmesher_extruder.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/** Creates nodes that are owned locally from the 2D template grid.*/
void chi_mesh::VolumeMesherExtruder::
CreateLocalNodes(chi_mesh::MeshContinuum& template_grid,
                 chi_mesh::MeshContinuum& grid)
{
  //================================================== For each layer
  std::set<uint64_t> vertex_ids_with_local_scope;
  for (size_t iz=0; iz<(vertex_layers.size()-1); iz++)
  {
    for (const auto& template_cell : template_grid.local_cells)
    {
      // Check template cell type
      if (template_cell.Type() != chi_mesh::CellType::POLYGON)
        throw std::logic_error("Extruder::CreateLocalNodes: "
                               "Template cell error. Not of base type POLYGON");

      bool has_local_scope = HasLocalScope(template_cell, template_grid, iz);

      if (has_local_scope)
      {
        auto& vertex_set = vertex_ids_with_local_scope;
        for (auto tc_vid : template_cell.vertex_ids)
          vertex_set.insert(tc_vid + iz * node_z_index_incr);

        for (auto tc_vid : template_cell.vertex_ids)
          vertex_set.insert(tc_vid + (iz + 1) * node_z_index_incr);
      }
    }//for template cell
  }//for layer

  //============================================= Now add all nodes
  //                                              that are local or neighboring
  uint64_t vid = 0;
  for (auto layer_z_level : vertex_layers)
  {
    for (auto& id_vertex : template_grid.vertices)
    {
      const auto& vertex = id_vertex.second;
      auto local_scope = vertex_ids_with_local_scope.find(vid);

      if (local_scope != vertex_ids_with_local_scope.end())
        grid.vertices.Insert(vid, Vector3(vertex.x, vertex.y, layer_z_level));

      ++vid;
    }//for vertex
  }//for layer

  grid.SetGlobalVertexCount(vid);
}