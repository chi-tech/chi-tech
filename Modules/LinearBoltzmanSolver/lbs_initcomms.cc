#include "lbs_linear_boltzman_solver.h"
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

//###################################################################
/**Initializes communicators*/
void LinearBoltzman::Solver::InitializeCommunicators()
{
  std::set<int>    local_graph_edges;
  std::vector<int> local_connections;

  //================================================== Loop over local cells
  //Populate local_graph_edges
  local_graph_edges.insert(chi_mpi.location_id); //add current location
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[c];
    auto cell           = grid->cells[cell_glob_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)cell;

      int num_faces = 2;
      for (int f=0; f<num_faces; f++)
      {
        int neighbor = slab_cell->edges[f];

        if (neighbor>=0)
        {
          auto adj_cell = grid->cells[neighbor];

          if (adj_cell->partition_id != chi_mpi.location_id)
          {
            local_graph_edges.insert(adj_cell->partition_id);
          }
        }
      }//for f
    } //if slab
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;

      for (int f=0; f< poly_cell->edges.size(); f++)
      {
        int neighbor = poly_cell->edges[f][EDGE_NEIGHBOR];

        if (neighbor>=0)
        {
          auto adj_cell = grid->cells[neighbor];

          if (adj_cell->partition_id != chi_mpi.location_id)
          {
            local_graph_edges.insert(adj_cell->partition_id);
          }
        }
      }//for f
    } //if polyhedron

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;

      for (int f=0; f< polyh_cell->faces.size(); f++)
      {
        int neighbor = polyh_cell->faces[f]->face_indices[NEIGHBOR];

        if (neighbor>=0)
        {
          auto adj_cell = grid->cells[neighbor];

          if (adj_cell->partition_id != chi_mpi.location_id)
          {
            local_graph_edges.insert(adj_cell->partition_id);
          }
        }
      }//for f
    } //if polyhedron
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encountered in call to "
           "InitializeCommunicators.";
      exit(EXIT_FAILURE);
    }
  }//for local cells

  //============================================= Convert set to vector
  //This is just done for convenience
  std::set<int>::iterator graph_edge;
  for (graph_edge =  local_graph_edges.begin();
       graph_edge != local_graph_edges.end();
       graph_edge++)
  {
    local_connections.push_back(*graph_edge);
  }

  //============================================= Broadcast local connection size
  chi_log.Log(LOG_0)
    << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(chi_mpi.process_count,
                                             std::vector<int>());
  for (int locI=0;locI<chi_mpi.process_count; locI++)
  {
    int locI_num_connections = local_connections.size();

    //If chi_mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(&locI_num_connections,1,MPI_INT,locI,MPI_COMM_WORLD);

    if (chi_mpi.location_id != locI)
    {global_graph[locI].resize(locI_num_connections,-1);}
    else
    {
      std::copy(local_connections.begin(),
                local_connections.end(),
                std::back_inserter(global_graph[locI]));
    }
  }



  //============================================= Broadcast local connections
  for (int locI=0;locI<chi_mpi.process_count; locI++)
  {
    //If chi_mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(global_graph[locI].data(),
              global_graph[locI].size(),
              MPI_INT,locI,MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_0)
    << "Done communicating local connections.";


  //============================================= Build groups
  MPI_Comm_group(MPI_COMM_WORLD,&comm_set.world_group);
  comm_set.location_groups.resize(chi_mpi.process_count,MPI_Group());

  for (int locI=0;locI<chi_mpi.process_count; locI++)
  {
    MPI_Group_incl(comm_set.world_group,
                   global_graph[locI].size(),
                   global_graph[locI].data(),
                   &comm_set.location_groups[locI]);
  }

  //============================================= Build communicators
  chi_log.Log(LOG_0)
    << "Building communicators.";
  comm_set.communicators.resize(chi_mpi.process_count,MPI_Comm());

  for (int locI=0;locI<chi_mpi.process_count; locI++)
  {
    int err = MPI_Comm_create_group(MPI_COMM_WORLD,
                                    comm_set.location_groups[locI],
                                    0, //tag
                                    &comm_set.communicators[locI]);

    if (!(err == MPI_SUCCESS))
    {
      chi_log.Log(LOG_ALL)
        << "Communicator creation failed.";
    }
  }

  chi_log.Log(LOG_0)
    << "Done building communicators.";

}