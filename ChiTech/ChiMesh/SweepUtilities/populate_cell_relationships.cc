#include "sweep_namespace.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_mpi.h>
#include <chi_log.h>


extern ChiLog& chi_log;

//###################################################################
/**Populates the local sub-grid connection information for sweep orderings.*/
void chi_mesh::sweep_management::PopulateCellRelationships(
         chi_mesh::MeshContinuumPtr grid,
         const chi_mesh::Vector3& omega,
         std::set<int>& location_dependencies,
         std::set<int>& location_successors,
         std::vector<std::set<std::pair<int,double>>>& cell_successors)
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
      double dot_normal = omega.Dot(face.normal);
      if (dot_normal>(0.0+tolerance)) {is_outgoing = true;}

      //======================================= If outgoing determine if
      //                                        it is to a local cell
      if (is_outgoing)
      {
        //================================ If it is a cell and not bndry
        if (face.has_neighbor)
        {
          //========================= If it is in the current location
          if (face.IsNeighborLocal(*grid))
          {
            double weight = dot_normal*face.ComputeFaceArea(*grid);
            cell_successors[c].insert(
              std::make_pair(face.GetNeighborLocalID(*grid),weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(*grid));
        }

      }
      //======================================= If not outgoing determine
      //                                        what it is dependent on
      else
      {
        //================================if it is a cell and not bndry
        if (face.has_neighbor and not face.IsNeighborLocal(*grid))
          location_dependencies.insert(face.GetNeighborPartitionID(*grid));
      }

    }//for edge
  }//for cell
}



