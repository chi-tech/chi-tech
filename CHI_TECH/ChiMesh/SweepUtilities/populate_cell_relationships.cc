#include "sweep_namespace.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

//###################################################################
/**Populates the local sub-grid connection information for sweep orderings.*/
void chi_mesh::sweep_management::PopulateCellRelationships(
         chi_mesh::MeshContinuum *grid,
         chi_mesh::sweep_management::SPDS* sweep_order,
         std::vector<std::set<int>>& cell_dependencies,
         std::vector<std::set<int>>& cell_successors)
{
  double tolerance = 1.0e-16;

  //============================================= Make directed connections
  for (auto& cell : grid->local_cells)
  {
    int c = cell.local_id;

    for (auto& face : cell.faces)
    {
      //======================================= Determine if the face
      //                                        is incident
      bool is_outgoing = false;
      double dot_normal = sweep_order->omega.Dot(face.normal);
      if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

      //======================================= If outgoing determine if
      //                                        it is to a local cell
      if (is_outgoing)
      {
        int adj_cell_glob_index = face.neighbor;

        //================================if it is a cell and not bndry
        if (adj_cell_glob_index>=0)
        {
          //========================= If it is in the current location
          if (face.IsNeighborLocal(grid))
          {
            cell_successors[c].insert(face.GetNeighborLocalID(grid));
          }
          else
            sweep_order->AddLocalSuccessor(face.GetNeighborPartitionID(grid));
        }

      }
        //======================================= If not outgoing determine
        //                                        what it is dependent on
      else
      {
        int adj_cell_glob_index = face.neighbor;

        //================================if it is a cell and not bndry
        if (adj_cell_glob_index>=0)
        {
          if (face.IsNeighborLocal(grid))
            cell_dependencies[c].insert(face.GetNeighborLocalID(grid));
          else
            sweep_order->AddLocalDependecy(face.GetNeighborPartitionID(grid));

        }
      }

    }//for edge
  }//for cell
}



