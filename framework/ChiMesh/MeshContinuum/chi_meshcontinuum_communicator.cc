#include "chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiMPI/chi_mpi_commset.h"

//###################################################################
/**Gets the communicator-set for interprocess communication,
 * associated with this mesh. If not created yet, it will create it.*/
chi_objects::ChiMPICommunicatorSet& chi_mesh::MeshContinuum::GetCommunicator()
{
  //================================================== Check if already avail
  if (communicator_set2_ != nullptr)
    return *communicator_set2_;

  //================================================== Build the communicator
  chi::log.Log0Verbose1() << "Building communicator.";
  std::set<int>    local_graph_edges;
  std::vector<int> local_connections;

  //================================================== Loop over local cells
  //Populate local_graph_edges
  local_graph_edges.insert(chi::mpi.location_id); //add current location
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor)
        if (not face.IsNeighborLocal(*this))
          local_graph_edges.insert(face.GetNeighborPartitionID(*this));
    }//for f
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
  chi::log.Log0Verbose1()
    << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(chi::mpi.process_count,
                                             std::vector<int>());
  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    int locI_num_connections = static_cast<int>(local_connections.size());

    //If chi::mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(&locI_num_connections,1,MPI_INT,locI,MPI_COMM_WORLD);

    if (chi::mpi.location_id != locI)
    {global_graph[locI].resize(locI_num_connections,-1);}
    else
    {
      std::copy(local_connections.begin(),
                local_connections.end(),
                std::back_inserter(global_graph[locI]));
    }
  }



  //============================================= Broadcast local connections
  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    //If chi::mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(global_graph[locI].data(),
              static_cast<int>(global_graph[locI].size()),
              MPI_INT,locI,MPI_COMM_WORLD);
  }

  chi::log.Log0Verbose1()
    << "Done communicating local connections.";


  //============================================= Build groups
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD,&world_group);

  std::vector<MPI_Group> location_groups;
  location_groups.resize(chi::mpi.process_count, MPI_Group());

  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    MPI_Group_incl(world_group,
                   static_cast<int>(global_graph[locI].size()),
                   global_graph[locI].data(),
                   &location_groups[locI]);
  }

  //============================================= Build communicators
  std::vector<MPI_Comm>  communicators;
  chi::log.Log0Verbose1()
    << "Building communicators.";
  communicators.resize(chi::mpi.process_count, MPI_Comm());

  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    int err = MPI_Comm_create_group(MPI_COMM_WORLD,
                                    location_groups[locI],
                                    0, //tag
                                    &communicators[locI]);

    if (err != MPI_SUCCESS)
    {
      chi::log.Log0Verbose1()
        << "Communicator creation failed.";
    }
  }

  chi::log.Log0Verbose1()
    << "Done building communicators.";

  communicator_set2_ = std::make_shared<chi_objects::ChiMPICommunicatorSet>(
    communicators, location_groups, world_group);

  return *communicator_set2_;
}

