#include "pwlc.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

#include <algorithm>

//###################################################################
/**Reorders the nodes for parallel computation in a Continuous
 * Finite Element calculation.*/
void SpatialDiscretization_PWLC::OrderNodes()
{
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Developing nodal ordering.";
  ChiTimer t_stage[6];

  t_stage[0].Reset();

  //================================================== Get SET of local
  //                                                   exclusive + non-exclusive
  //                                                   nodes
  std::set<int> exnonex_nodes_set;
  for (auto& cell : ref_grid->local_cells)
    for (auto& vid : cell.vertex_ids)
      exnonex_nodes_set.insert(vid);

  // Copy set into vector
  std::vector<int> exnonex_nodes;
  for (auto& vid : exnonex_nodes_set)
    exnonex_nodes.push_back(vid);

  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 0 time: "
                              << t_stage[0].GetTime()/1000.0;

  t_stage[1].Reset();
  // Initialize ghost flags
  std::vector<bool> ghost_flags;
  ghost_flags.resize(exnonex_nodes.size(),false);


  //================================================== Determine ghost flags
  //Run through each cell and for each local nodes deterime
  //if the nodes are ghosts to another location
  for (auto& cell : ref_grid->local_cells)
  {
    for (size_t f=0; f < cell.faces.size(); f++)
    {
      auto& face = cell.faces[f];
      if (face.has_neighbor)
      {
        if (face.GetNeighborPartitionID(*ref_grid) != cell.partition_id)
        {
          for (auto v_index : face.vertex_ids)
          {
            int v_setind =
              (int)std::distance(exnonex_nodes.begin(),
                                 std::find(std::begin(exnonex_nodes),
                                           std::end(exnonex_nodes),
                                           v_index));
            ghost_flags[v_setind] = true;
          }//for face verts
        }//if neighbor not local
      }//if neighbor is not a boundary
    }//for cell face
  }
  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 1 time: "
                              << t_stage[1].GetTime()/1000.0;

  t_stage[2].Reset();
  //================================================== Develop vectors of exclusive
  //                                                   and non-exclusive nodes
  // Basically splitting exnonex_nodes into exclusive
  // and non-exclusive nodes
  std::vector<int> exclusive_nodes;
  std::vector<int> nonexclus_nodes;
  for (size_t i=0; i<exnonex_nodes.size(); i++)
  {
    int ind = exnonex_nodes[i];
    if (ghost_flags[i])
      nonexclus_nodes.push_back(ind);
    else
      exclusive_nodes.push_back(ind);
  }

  chi_log.Log(LOG_ALLVERBOSE_1) << "Number of exclusive nodes = "
                                << exclusive_nodes.size()
                                << " and ghost nodes = "
                                << nonexclus_nodes.size();
  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 2 time: "
                              << t_stage[2].GetTime()/1000.0;

  t_stage[3].Reset();
  //================================================== Ring communicate ghost nodes
  //Using ring communication we communicate the ghost
  //nodes to find a single list of all the ghost nodes
  std::vector<int> global_ghost_nodes;

  // Location 0 sends to location 1
  if (chi_mpi.location_id==0)
  {
    MPI_Send(nonexclus_nodes.data(),
             nonexclus_nodes.size(),
             MPI_INT,1,123,MPI_COMM_WORLD);
  }
    // Location n=1..(N-1) first receives n-1
  else
  {
    std::vector<int> upstream_nonex;

    MPI_Status status;
    MPI_Probe(chi_mpi.location_id-1,123,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    upstream_nonex.resize(num_to_recv,-1);
    MPI_Recv(upstream_nonex.data(),num_to_recv,
             MPI_INT,chi_mpi.location_id-1,123,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Run through location n non-exclusive nodes
    // if a non-exclusive node is not in the upstream list then
    // it gets added to the list
    //for (size_t i=0; i<nonexclus_nodes.size(); i++)
    for (auto& node : nonexclus_nodes)
    {
      auto node_location = std::find(upstream_nonex.begin(),
                                     upstream_nonex.end(),
                                     node);
      if (node_location == upstream_nonex.end())
        upstream_nonex.push_back(node);
    }

    // If this location is not the last location
    // send the updated upstream_nonex to location n+1
    if (chi_mpi.location_id<(chi_mpi.process_count-1))
      MPI_Send(upstream_nonex.data(),upstream_nonex.size(),
               MPI_INT,chi_mpi.location_id+1,123,MPI_COMM_WORLD);
      // On the last location send the completed
      // upstream_nonex back to all other locations
    else
    {
      chi_log.Log(LOG_ALLVERBOSE_1)
        << "Total number of ghost nodes after collect: "
        << upstream_nonex.size();
      std::copy(upstream_nonex.begin(),
                upstream_nonex.end(),
                std::back_inserter(global_ghost_nodes));

      for (int loc=0; loc<(chi_mpi.process_count-1); loc++)
        MPI_Send(global_ghost_nodes.data(),global_ghost_nodes.size(),
                 MPI_INT,loc,124,MPI_COMM_WORLD);

    }
  }

  //================================================== Last location broadcasts
  if (chi_mpi.location_id<(chi_mpi.process_count-1))
  {
    MPI_Status status;
    MPI_Probe(chi_mpi.process_count-1,124,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    global_ghost_nodes.resize(num_to_recv,-1);
    MPI_Recv(global_ghost_nodes.data(),num_to_recv,
             MPI_INT,chi_mpi.process_count-1,124,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  chi_log.Log(LOG_ALLVERBOSE_1) << "Total number of ghost nodes: "
                                << global_ghost_nodes.size() << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 3 time: "
                              << t_stage[3].GetTime()/1000.0;

  t_stage[3].Reset();


  //================================================== Seperate ghost nodes into
  //                                                   pieces
  // Imagine that each location is a cell in a 1D array.
  // The interfaces between each cell is where the ghost nodes are
  // and therefore each cell needs to take ownership of the some ghost nodes.
  // Now the ghost nodes to the left of an interior cell is the same
  // as the ghost nodes to the right of cell i-1. One approach is then
  // to split these interface ghost nodes in half and give each cell
  // half of them. Therefore we split each interface in half and assign
  // the ghost node-halfs to cells on either side of the interface.
  // With N amount of "cells" there would be N-1 interfaces and therefore
  // 2(N-1) halfs of the ghost nodes.
  int g_piece = 0;
  if (chi_mpi.process_count>1)
    g_piece = (int)floor(global_ghost_nodes.size()/
                         (2*(chi_mpi.process_count-1)));

  int g_from = 0;
  int g_to = 0;
  if (chi_mpi.location_id==0)
  {
    g_from = 0;
    g_to = g_piece-1;
  }
  else if (chi_mpi.location_id < (chi_mpi.process_count-1))
  {
    g_from = g_piece-1 + 2*g_piece*(chi_mpi.location_id-1)+1;
    g_to   = g_from + 2*g_piece-1;
  }
  else
  {
    g_from = g_piece-1 + 2*g_piece*(chi_mpi.location_id-1)+1;
    g_to   = (int)global_ghost_nodes.size()-1;
  }
  int num_g_loc = (g_to - g_from +1);
  chi_log.Log(LOG_ALLVERBOSE_1) << "Local ghost ownership: "
                                << g_from << "->" << g_to
                                << "(" << num_g_loc << ")" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Ring Compute local portion
  //The local portion of the nodes are the exclusive nodes
  //plus the locally owned ghost nodes
  //The list of exclusive nodes will be packed into
  //a parallel array according it ownership. Lets say
  //each location n has ownership of the parallel vectors starting
  //at x_ni. x_0i should obviously be x_0i=0. Now, location n has ownership
  //from x_ni to x_nf where we still have to determine x_f. Somewhere between
  //x_ni and x_nf there would be x_n_local (the point where the exclusive nodes
  //ends and where we need to pack the portion of the ghost nodes that location
  // n took ownership of as determined by the previous step.
  int local_from = 0;
  int local_to = 0;
  if (chi_mpi.location_id==0)
  {
    local_from = 0;
    local_to = (int)exclusive_nodes.size() - 1 + num_g_loc;
    MPI_Send(&local_to,1,
             MPI_INT,chi_mpi.location_id+1,125,MPI_COMM_WORLD);
  }
  else
  {
    int upstream_loc_end = 0;
    MPI_Recv(&upstream_loc_end,1,
             MPI_INT,chi_mpi.location_id-1,125,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    local_from = upstream_loc_end + 1;
    local_to = local_from + (int)exclusive_nodes.size() - 1 + num_g_loc;

    if (chi_mpi.location_id<(chi_mpi.process_count-1))
      MPI_Send(&local_to,1,
               MPI_INT,chi_mpi.location_id+1,125,MPI_COMM_WORLD);

  }
  int tot_local_nodes = local_to - local_from + 1;
  chi_log.Log(LOG_ALLVERBOSE_1) << "Local node ownership "
                                << local_from << "->" << local_to
                                << "(" << tot_local_nodes << ")"
                                << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 4 time: "
                              << t_stage[4].GetTime()/1000.0;
  MPI_Barrier(MPI_COMM_WORLD);

  t_stage[5].Reset();


  //================================================== Ring map the ghost nodes
  // We still have a list of global_ghost_nodes that needs mapping but
  // fortunately we know which portions of this we own [g_from to g_to].
  // We can now compute a mapping of the ghost nodes to where they will be
  // placed in the global array. The mapping will not be complete yet
  // since it is still based on the ghost node index in the vector
  // global_ghost_nodes.
  std::vector<int> ghost_mapping;
  ghost_mapping.resize(global_ghost_nodes.size(),-1);
  if (chi_mpi.location_id==0)
  {
    for (int g=g_from; g<= g_to; g++)
    {
      int ghost_index = local_from + (int)exclusive_nodes.size() + g-g_from;
      ghost_mapping[g] = ghost_index;
    }
    MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
             MPI_INT,chi_mpi.location_id+1,126,MPI_COMM_WORLD);
  }
  else
  {
    std::vector<int> upstream_ghost_mapping;
    MPI_Status status;
    MPI_Probe(chi_mpi.location_id-1,126,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    upstream_ghost_mapping.resize(num_to_recv,-1);
    MPI_Recv(upstream_ghost_mapping.data(),num_to_recv,
             MPI_INT,chi_mpi.location_id-1,126,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    std::copy(upstream_ghost_mapping.begin(),
              upstream_ghost_mapping.end(),
              ghost_mapping.begin());
    for (int g=g_from; g<= g_to; g++)
    {
      int ghost_index = local_from + (int)exclusive_nodes.size() + g-g_from;
      ghost_mapping[g] = ghost_index;
    }

    if (chi_mpi.location_id<(chi_mpi.process_count-1))
      MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
               MPI_INT,chi_mpi.location_id+1,126,MPI_COMM_WORLD);
    else
      for (int loc=0; loc<(chi_mpi.process_count-1); loc++)
        MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
                 MPI_INT,loc,127,MPI_COMM_WORLD);
  }

  //================================================== Collect ghost mapping
  //                                                   from last location
  if (chi_mpi.location_id<(chi_mpi.process_count-1))
  {
    MPI_Status status;
    MPI_Probe(chi_mpi.process_count-1,127,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    ghost_mapping.resize(num_to_recv,-1);
    MPI_Recv(ghost_mapping.data(),num_to_recv,
             MPI_INT,chi_mpi.process_count-1,127,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }


  //================================================== Initialize forward
  //                                                   ordering
  node_mapping.clear();


  //================================================== Creating mapping of
  //                                                   exclusive nodes
  int num_exl = exclusive_nodes.size();
  for (int n=0; n<num_exl; n++)
  {
    int orig_index = exclusive_nodes[n];
    node_mapping[orig_index] = local_from + n;
  }


  //================================================== Create mapping of ghost
  //                                                   nodes
  for (size_t n=0; n<global_ghost_nodes.size(); n++)
  {
    int orig_index = global_ghost_nodes[n];
    node_mapping[orig_index] = ghost_mapping[n];
  }

  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 5 time: "
                              << t_stage[5].GetTime()/1000.0;
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Compute block addresses
  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stages complete time: "
                              << t_stage[5].GetTime()/1000.0;
  MPI_Barrier(MPI_COMM_WORLD);

  local_block_address = local_from;

  local_base_block_size = local_to - local_from + 1;
  globl_base_block_size = ref_grid->GetGlobalVertexCount();

  //======================================== Collect block addresses
  locJ_block_address.clear();
  locJ_block_address.resize(chi_mpi.process_count, 0);
  MPI_Allgather(&local_block_address,    //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_block_address.data(),    //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator

  //======================================== Collect block sizes
  locJ_block_size.clear();
  locJ_block_size.resize(chi_mpi.process_count, 0);
  MPI_Allgather(&local_base_block_size,       //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_block_size.data(),       //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator
}