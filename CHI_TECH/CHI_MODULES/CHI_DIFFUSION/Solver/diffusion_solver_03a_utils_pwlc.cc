#include "diffusion_solver.h"
#include "../../../ChiMesh/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../ChiMesh/CHI_REGION/chi_region.h"
#include "../../../ChiMesh/CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../../../ChiMesh/CHI_CELL/cell.h"
#include "../../../ChiMesh/CHI_CELL/cell_slab.h"
#include "../../../ChiMesh/CHI_CELL/cell_polygon.h"
#include "../../../ChiMesh/CHI_CELL/cell_polyhedron.h"
#include "../../../ChiMesh/CHI_VOLUMEMESHER/chi_volumemesher.h"
#include "../../../ChiTimer/chi_timer.h"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;



//###################################################################
/**Reorders nodes for better parrallel communication during matrix
 * assembly.
 *
 * The first stage of reordering needs to make an attempt at ordering
 * nodes that are mostly local to the local cells.*/
void chi_diffusion::Solver::ReorderNodesPWLC()
{
  ChiTimer t_stage[6];

  t_stage[0].Reset();
  //================================================== Get reference to continuum
  auto handler = chi_mesh::GetCurrentHandler();
  auto region  = handler->region_stack.back();
  auto vol_continuum = region->volume_mesh_continua.back();

  auto mesher = handler->volume_mesher;

  //================================================== Get SET of local
  //                                                   exclusive + non-exclusive
  //                                                   nodes
  std::set<int> exnonex_nodes_set;
  size_t num_loc_cells = vol_continuum->local_cell_glob_indices.size();
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = vol_continuum->local_cell_glob_indices[lc];
    auto cell = vol_continuum->cells[cell_glob_index];

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      for (int v=0; v<2; v++)
      {
        exnonex_nodes_set.insert(slab_cell->v_indices[v]);
      }
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      for (int v=0; v<poly_cell->v_indices.size(); v++)
      {
        exnonex_nodes_set.insert(poly_cell->v_indices[v]);
      }
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      for (int v=0; v<polyh_cell->v_indices.size(); v++)
      {
        exnonex_nodes_set.insert(polyh_cell->v_indices[v]);
      }
    }
  }

  //================================================== Copy set into vector
  std::vector<int> exnonex_nodes;
  std::set<int>::iterator set_iter;
  for (set_iter=exnonex_nodes_set.begin();
       set_iter != exnonex_nodes_set.end(); set_iter++)
  {
    exnonex_nodes.push_back(*set_iter);
  }
  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 0 time: "
                            << t_stage[0].GetTime()/1000.0;

  t_stage[1].Reset();
  //================================================== Initialize ghost flags
  std::vector<bool> ghost_flags;
  ghost_flags.resize(exnonex_nodes.size(),false);


  //================================================== Determine ghost flags
  //Run through each cell and for each local nodes deterime
  //if the nodes are ghosts to another location
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = vol_continuum->local_cell_glob_indices[lc];
    auto cell = vol_continuum->cells[cell_glob_index];

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      for (int e=0; e<2; e++)
      {
        if (slab_cell->edges[e]>=0)
        {
          int adj_cell_ind = slab_cell->edges[e];
          auto adj_cell = vol_continuum->cells[adj_cell_ind];

          if (adj_cell->partition_id != slab_cell->partition_id)
          {
            for (int ev=0; ev<2; ev++)
            {
              int v_index = slab_cell->edges[ev];
              int v_setind = (int)std::distance(exnonex_nodes.begin(),
                                                std::find(std::begin(exnonex_nodes),
                                                          std::end(exnonex_nodes),
                                                          v_index));
              ghost_flags[v_setind] = true;
            }//for edge vertices
          }//if neighbor not local
        }//if neigbor not a boundary
      }//for each edge
    }//if slab

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      for (int e=0; e<poly_cell->edges.size(); e++)
      {
        if (poly_cell->edges[e][2]>=0)
        {
          int adj_cell_ind = poly_cell->edges[e][2];
          auto adj_cell = vol_continuum->cells[adj_cell_ind];

          if (adj_cell->partition_id != poly_cell->partition_id)
          {
            for (int ev=0; ev<2; ev++)
            {
              int v_index = poly_cell->edges[e][ev];
              int v_setind = (int)std::distance(exnonex_nodes.begin(),
                                                std::find(std::begin(exnonex_nodes),
                                                          std::end(exnonex_nodes),
                                                          v_index));
              ghost_flags[v_setind] = true;
            }//for edge vertices
          }//if neighbor not local
        }//if neigbor not a boundary
      }//for each edge
    }//if polygon

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      for (int f=0; f<polyh_cell->faces.size(); f++)
      {
        if (polyh_cell->faces[f]->face_indices[NEIGHBOR]>=0)
        {
          int adj_cell_ind = polyh_cell->faces[f]->face_indices[0];
          auto adj_cell = vol_continuum->cells[adj_cell_ind];

          if (adj_cell->partition_id != polyh_cell->partition_id)
          {
            for (int fv=0; fv<polyh_cell->faces[f]->v_indices.size(); fv++)
            {
              int v_index = polyh_cell->faces[f]->v_indices[fv];
              int v_setind = (int)std::distance(exnonex_nodes.begin(),
                                           std::find(std::begin(exnonex_nodes),
                                                     std::end(exnonex_nodes),
                                                     v_index));
              ghost_flags[v_setind] = true;
            }//for face verts
          }//if neighbor not local
        }//if neighbor is not a boundary
      }//for cell face
    }//if polyhedron
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
  for (int i=0; i<exnonex_nodes.size(); i++)
  {
    int ind = exnonex_nodes[i];//*std::next(exnonex_nodes.begin(),i);
    if (ghost_flags[i])
      nonexclus_nodes.push_back(ind);
    else
      exclusive_nodes.push_back(ind);

    //chi_log.Log(LOG_0) << i << " " << ind;
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
//  mpi::environment env;
//  mpi::communicator world;

  std::vector<int> global_ghost_nodes;

  //=================================== Location 0 sends to location 1
  if (chi_mpi.location_id==0)
  {
    //world.send(1,123,nonexclus_nodes);
    MPI_Send(nonexclus_nodes.data(),nonexclus_nodes.size(),
             MPI_INT,1,123,MPI_COMM_WORLD);
  }
  //=================================== Location n=1..(N-1) first receives
  //                                    n-1
  else
  {
    std::vector<int> upstream_nonex;
    //world.recv(chi_mpi.location_id-1,123,upstream_nonex);
    MPI_Status status;
    MPI_Probe(chi_mpi.location_id-1,123,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    upstream_nonex.resize(num_to_recv,-1);
    MPI_Recv(upstream_nonex.data(),num_to_recv,
             MPI_INT,chi_mpi.location_id-1,123,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    //============================ Run through location n non-exclusive nodes
    // if a non-exclusive node is not in the upstream list then
    // it gets added to the list
    for (int i=0; i<nonexclus_nodes.size(); i++)
    {
      if (std::find(upstream_nonex.begin(),
                    upstream_nonex.end(),
                    nonexclus_nodes[i]) == upstream_nonex.end())
        upstream_nonex.push_back(nonexclus_nodes[i]);
    }

    //============================ If this location is not the last location
    // send the updated upstream_nonex to location n+1
    if (chi_mpi.location_id<(chi_mpi.process_count-1))
    {
      //world.send(chi_mpi.location_id+1,123,upstream_nonex);
      MPI_Send(upstream_nonex.data(),upstream_nonex.size(),
               MPI_INT,chi_mpi.location_id+1,123,MPI_COMM_WORLD);
    }
    //============================ On the last location send the completed
    // upstream_nonex back to all other locations
    else
    {
      chi_log.Log(LOG_ALLVERBOSE_1) << "Total number of ghost nodes after collect: "
                                << upstream_nonex.size();
      std::copy(upstream_nonex.begin(),upstream_nonex.end(),std::back_inserter(global_ghost_nodes));
      for (int loc=0; loc<(chi_mpi.process_count-1); loc++)
      {
//        world.send(loc,124,global_ghost_nodes);
        MPI_Send(global_ghost_nodes.data(),global_ghost_nodes.size(),
                 MPI_INT,loc,124,MPI_COMM_WORLD);
      }

    }
  }



  //================================================== Collect ghost nodes
  if (chi_mpi.location_id<(chi_mpi.process_count-1))
  {
//    world.recv(chi_mpi.process_count-1,124,global_ghost_nodes);
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
  {
    g_piece = (int)floor(global_ghost_nodes.size()/
                         (2*(chi_mpi.process_count-1)));
  }

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
  } else
  {
    g_from = g_piece-1 + 2*g_piece*(chi_mpi.location_id-1)+1;
    g_to   = global_ghost_nodes.size()-1;
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
    local_to = exclusive_nodes.size() - 1 + num_g_loc;
    //world.send(chi_mpi.location_id+1,125,local_to);
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
    local_to = local_from + exclusive_nodes.size() - 1 + num_g_loc;

    if (chi_mpi.location_id<(chi_mpi.process_count-1))
    {
      MPI_Send(&local_to,1,
               MPI_INT,chi_mpi.location_id+1,125,MPI_COMM_WORLD);
    }

  }
  int tot_local_nodes = local_to - local_from + 1;
  chi_log.Log(LOG_ALLVERBOSE_1) << "Local node ownership "
                                      << local_from << "->" << local_to
                                      << "(" << tot_local_nodes << ")"
                                      << std::endl;
//  usleep(1000000);
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
      int ghost_index = local_from + exclusive_nodes.size() + g-g_from;
      ghost_mapping[g] = ghost_index;
    }
    //world.send(chi_mpi.location_id+1,126,ghost_mapping);
    MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
             MPI_INT,chi_mpi.location_id+1,126,MPI_COMM_WORLD);
  }
  else
  {
    std::vector<int> upstream_ghost_mapping;
    //world.recv(chi_mpi.location_id-1,126,upstream_ghost_mapping);
    MPI_Status status;
    MPI_Probe(chi_mpi.location_id-1,126,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    upstream_ghost_mapping.resize(num_to_recv,-1);
    MPI_Recv(upstream_ghost_mapping.data(),num_to_recv,
             MPI_INT,chi_mpi.location_id-1,126,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);

//    ghost_mapping = upstream_ghost_mapping;
    std::copy(upstream_ghost_mapping.begin(),
              upstream_ghost_mapping.end(),
              ghost_mapping.begin());
    for (int g=g_from; g<= g_to; g++)
    {
      int ghost_index = local_from + exclusive_nodes.size() + g-g_from;
      ghost_mapping[g] = ghost_index;
    }
    if (chi_mpi.location_id<(chi_mpi.process_count-1))
    {
//      world.send(chi_mpi.location_id+1,126,ghost_mapping);
      MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
               MPI_INT,chi_mpi.location_id+1,126,MPI_COMM_WORLD);
    }
    else
    {
      for (int loc=0; loc<(chi_mpi.process_count-1); loc++)
      {
        //world.send(loc,127,ghost_mapping);//
        MPI_Send(ghost_mapping.data(),ghost_mapping.size(),
                 MPI_INT,loc,127,MPI_COMM_WORLD);
      }
    }
  }

  //================================================== Collect ghost mapping
  //                                                   from last location
  if (chi_mpi.location_id<(chi_mpi.process_count-1))
  {
//    world.recv(chi_mpi.process_count-1,127,ghost_mapping);
    MPI_Status status;
    MPI_Probe(chi_mpi.process_count-1,127,MPI_COMM_WORLD,&status);
    int num_to_recv=0;
    MPI_Get_count(&status,MPI_INT,&num_to_recv);
    ghost_mapping.resize(num_to_recv,-1);
    MPI_Recv(ghost_mapping.data(),num_to_recv,
             MPI_INT,chi_mpi.process_count-1,127,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }


  //$$$$$$$$$$$$$$$ TEMP CODE
//  if (chi_mpi.location_id==2)
//  {
//    FILE* ofile = std::fopen("Loc2Ghost.csv","w");
//    printf("Writing Loc2\n");
//    for (int g=0; g<ghost_mapping.size(); g++)
//    {
//      fprintf(ofile,"%d %d\n",g,ghost_mapping[g]);
//    }
//    std::fclose(ofile);
//  }
  //$$$$$$$$$$$$$$$ END TEMP CODE


  //================================================== Initialize forward
  //                                                   ordering
  //Initially this is just pass through
  std::vector<int> node_ordering;
  int num_nodes = vol_continuum->nodes.size();
  node_ordering.resize(num_nodes,-1);


  //================================================== Creating mapping of
  //                                                   exclusive nodes
  int num_exl = exclusive_nodes.size();
  for (int n=0; n<num_exl; n++)
  {
    int orig_index = exclusive_nodes[n];
    node_ordering[orig_index] = local_from + n;
  }


  //================================================== Create mapping of ghost
  //                                                   nodes
  for (int n=0; n<global_ghost_nodes.size(); n++)
  {
    int orig_index = global_ghost_nodes[n];
    node_ordering[orig_index] = ghost_mapping[n];
  }

  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stage 5 time: "
                            << t_stage[5].GetTime()/1000.0;
  MPI_Barrier(MPI_COMM_WORLD);


  //================================================== Push up these mappings
  //                                                   to the mesher
  mesher->node_ordering.clear();
  mesher->reverse_node_ordering.clear();
  mesher->reverse_node_ordering.resize(num_nodes,-1);
  int temp;
  for (int i=0; i<num_nodes; i++)
  {
    mesher->node_ordering.push_back(new chi_mesh::NodeIndexMap(i,node_ordering[i]));
    temp = node_ordering[i];

    if (temp>=0)
    {
      mesher->reverse_node_ordering[node_ordering[i]] = i;
    }
  }

  this->local_rows_from = local_from;
  this->local_rows_to = local_to;

  chi_log.Log(LOG_0VERBOSE_1) << "*** Reordering stages complete time: "
                            << t_stage[5].GetTime()/1000.0;
  MPI_Barrier(MPI_COMM_WORLD);
//  exit(0);
}


