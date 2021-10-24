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
  std::set<int> local_vert_ids;
  for (int iz=0; iz<(vertex_layers.size()-1); iz++)
  {
    for (int tc=0; tc < template_grid.local_cells.size(); tc++)
    {
      //========================================= Get template cell
      if (template_grid.local_cells[tc].Type() !=
          chi_mesh::CellType::POLYGON)
      {
        chi_log.Log(LOG_ALLERROR)
          << "Extruder::CreateLocalAndBoundaryNodes: Template cell error.";
        exit(EXIT_FAILURE);
      }
      auto& template_cell = template_grid.local_cells[tc];

      //========================================= Precompute centroid
      auto centroid_precompd = ComputeTemplateCell3DCentroid(
        template_cell, template_grid, iz, iz + 1);

      //========================================= Get the partition id
      int tcell_partition_id = GetCellPartitionIDFromCentroid(centroid_precompd);

      //###################### NOT A LOCAL CELL ############################
      if (tcell_partition_id != chi_mpi.location_id)
      {
        bool is_neighbor_to_partition = IsTemplateCellNeighborToThisPartition(
          template_cell, template_grid, iz, tc);

        if (is_neighbor_to_partition)
        {
          for (auto tc_vid : template_cell.vertex_ids)
            local_vert_ids.insert(tc_vid + iz*node_z_index_incr);

          for (auto tc_vid : template_cell.vertex_ids)
            local_vert_ids.insert(tc_vid + (iz+1)*node_z_index_incr);
        }
      }
        //####################### LOCAL CELL ###############################
      else
      {
        for (auto tc_vid : template_cell.vertex_ids)
          local_vert_ids.insert(tc_vid + iz*node_z_index_incr);

        for (auto tc_vid : template_cell.vertex_ids)
          local_vert_ids.insert(tc_vid + (iz+1)*node_z_index_incr);
      }
    }//for template cell
  }//for layer

  //============================================= Now add all nodes
  //                                              that are local or neighboring
  grid.vertices.clear();

  for (auto layer_z_level : vertex_layers)
  {
    for (auto& vertex : template_grid.vertices)
    {
      size_t new_vert_index = grid.vertices.size();

      auto local_index = local_vert_ids.find(new_vert_index);

      if (local_index != local_vert_ids.end())
      {
        auto node = vertex;
        node.z = layer_z_level;

        grid.vertices.push_back(node);
      }
      else
        grid.vertices.emplace_back();

    }//for vertex
  }//for layer
}