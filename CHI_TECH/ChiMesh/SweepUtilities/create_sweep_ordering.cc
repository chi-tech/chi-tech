#include "../chi_mesh.h"
#include "sweep_namespace.h"

#include "../MeshHandler/chi_meshhandler.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include "../../ChiConsole/chi_console.h"
#include <ChiTimer/chi_timer.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;
extern ChiConsole chi_console;
extern ChiTimer   chi_program_timer;

#include <ChiGraph/chi_directed_graph.h>

//###################################################################
/**Develops a sweep ordering for a given angle for locally owned
 * cells.*/
chi_mesh::sweep_management::SPDS* chi_mesh::sweep_management::
  CreateSweepOrder(double polar, double azimuthal,
                   chi_mesh::MeshContinuum *grid,
                   bool cycle_allowance_flag)
{
  auto sweep_order  = new chi_mesh::sweep_management::SPDS;
  sweep_order->grid = grid;

  size_t num_loc_cells = grid->local_cell_glob_indices.size();

  //============================================= Compute direction vector
  sweep_order->polar     = polar;
  sweep_order->azimuthal = azimuthal;

  sweep_order->omega.x = sin(polar)*cos(azimuthal);
  sweep_order->omega.y = sin(polar)*sin(azimuthal);
  sweep_order->omega.z = cos(polar);

  chi_mesh::Vector3 omega = sweep_order->omega; //shorter name
  if (chi_mpi.location_id == 0)
  {
    char buff[100];
    snprintf(buff, sizeof(buff), "Omega = %f,%f,%f\n",omega.x,omega.y,omega.z);

    chi_log.Log(LOG_0VERBOSE_1) << buff;
  }

  //============================================= Make directed connections
  chi_log.Log(LOG_0VERBOSE_1) << "Populating cell relationships";
  std::vector<std::set<int>> cell_dependencies(num_loc_cells);
  std::vector<std::set<int>> cell_successors(num_loc_cells);
  PopulateCellRelationships(grid,
                            sweep_order,
                            cell_dependencies,
                            cell_successors);

  //============================================= Build graph
  chi_graph::DirectedGraph test_DG;

  for (int c=0; c<num_loc_cells; ++c)
    test_DG.AddVertex();

  for (int c=0; c<num_loc_cells; c++)
    for (auto successor : cell_successors[c])
      test_DG.AddEdge(c,successor);

  //============================================= Find initial SCCs
  auto SCCs = test_DG.FindStronglyConnectedConnectionns();

  //============================================= Rinse remove edges
  std::vector<std::pair<int,int>> edges_to_remove;
  int iter=0;
  while (not SCCs.empty())
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "Local cyclic dependency removal. Iteration " << ++iter;
    //=================================== Loop over sub-graphs
    edges_to_remove.clear();
    for (auto& subDG : SCCs)
    {
      //Identify edges to remove
      bool found_edge_to_remove = false;
      for (int v : subDG)
      {
        for (int w : test_DG.vertices[v].ds_edge)
          if (std::find(subDG.begin(), subDG.end(), w) != subDG.end())
          {
            found_edge_to_remove = true;
            edges_to_remove.emplace_back(v,w);
            break;
          }

        if (found_edge_to_remove) break;
      }//for v in subDG
    }//for subDG

    //Remove the edges
    for (auto& edge_to_remove : edges_to_remove)
    {
      test_DG.RemoveEdge(edge_to_remove.first,edge_to_remove.second);
      sweep_order->local_cyclic_dependencies.emplace_back(
                                               edge_to_remove.first,
                                               edge_to_remove.second);
    }

    // Refind SCCs
    SCCs = test_DG.FindStronglyConnectedConnectionns();
  }

  //============================================= Generate topological sorting
  chi_log.Log(LOG_0VERBOSE_1) << "Generating topological sorting";
  sweep_order->spls = new chi_mesh::sweep_management::SPLS;
  sweep_order->spls->item_id = test_DG.GenerateTopologicalSort();

  if (sweep_order->spls->item_id.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Topological sorting for local sweep-ordering failed. "
      << "Cyclic dependencies detected.";
    exit(EXIT_FAILURE);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Task Dependency Graphs
  //All locations will send their dependencies
  //to the other locations

  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Communicating sweep dependencies.";

  int P = chi_mpi.process_count;
  std::vector<int> dependency_count_per_location(P,0);
  std::vector<std::vector<int>> global_dependencies(P, std::vector<int>());

  dependency_count_per_location[chi_mpi.location_id] =
    sweep_order->location_dependencies.size();

  //============================================= Broadcast location dep counts
  for (int locI=0; locI<P; locI++)
    MPI_Bcast(&dependency_count_per_location[locI], //Buffer
              1, MPI_INT,                           //Count and type
              locI,                                 //Sending location
              MPI_COMM_WORLD);                      //Communicator

  //============================================= Broadcast dependencies
  for (int locI=0; locI<P; locI++)
  {
    if (locI == chi_mpi.location_id)
    {
      std::copy(sweep_order->location_dependencies.begin(),
                sweep_order->location_dependencies.end(),
                std::back_inserter(global_dependencies[locI]));
    }
    else
    {
      global_dependencies[locI].
        resize(dependency_count_per_location[locI],-1);
    }

    MPI_Bcast(global_dependencies[locI].data(),        //Buffer
              dependency_count_per_location[locI],     //Count
              MPI_INT,                                 //Type
              locI,                                    //Sending location
              MPI_COMM_WORLD);                         //Communicator
  }

  //====================================== Filter dependencies for cycles
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Removing intra-cellset cycles.";

  const bool ALLOW_RECURSIVE_SEARCH = true;
  const bool DONT_ALLOW_RECURSIVE_SEARCH = false;

  RemoveGlobalCyclicDependencies(sweep_order, global_dependencies,
                                 DONT_ALLOW_RECURSIVE_SEARCH, cycle_allowance_flag);
  RemoveGlobalCyclicDependencies(sweep_order, global_dependencies,
                                 ALLOW_RECURSIVE_SEARCH, cycle_allowance_flag);
  MPI_Barrier(MPI_COMM_WORLD);

  //====================================== Build task dependency graph
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Building Task Dependency Graphs.";
  CHI_D_GRAPH TDG;
  std::vector<int> glob_linear_sweep_order;
  std::vector<int> glob_sweep_order_rank;
  std::vector<int> glob_order_mapping;

  //================================= Add vertices to the graph
  for (int loc=0; loc<chi_mpi.process_count; loc++)
  {
    boost::add_vertex(TDG);
    glob_sweep_order_rank.push_back(-1);
    glob_order_mapping.push_back(loc);
  }

  chi_log.Log(LOG_ALLVERBOSE_1)
    << chi_program_timer.GetTimeString()
    << "   - Adding dependencies.";
  //================================= Add dependencies
  for (int loc=0; loc<chi_mpi.process_count; loc++)
  {
    chi_log.Log(LOG_0VERBOSE_1) << "Location " << loc;
    for (int dep=0; dep<global_dependencies[loc].size(); dep++)
    {
      if (global_dependencies[loc][dep]>=0)
      {
        boost::add_edge(global_dependencies[loc][dep], loc, TDG);

        chi_log.Log(LOG_0VERBOSE_1) << " depends loc " << global_dependencies[loc][dep];
      }

    }
  }

  //================================= Generate topological sort
  chi_log.Log(LOG_ALLVERBOSE_1)
    << chi_program_timer.GetTimeString()
    << "   - Generating topological sort.";
  typedef boost::graph_traits<CHI_D_GRAPH>::vertex_descriptor gVertex;

  boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
    glob_index_map = get(boost::vertex_index, TDG);

  std::vector<gVertex> glob_sorted_list;
  try {
    boost::topological_sort(TDG, std::back_inserter(glob_sorted_list));
  }
  catch (const boost::bad_graph& exc)
//  catch(...)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Cyclic global sweep ordering detected.";
    exit(EXIT_FAILURE);
  }

  //================================= Generate linear ordering
  // This step just puts the topological
  // sorting into a std::vector.
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Generating global ordering.";
  for (std::vector<gVertex>::reverse_iterator ii=glob_sorted_list.rbegin();
       ii!=glob_sorted_list.rend();
       ii++)
  {
    int location_index = glob_index_map[*ii];
    glob_linear_sweep_order.push_back(location_index);
  }
  if (glob_linear_sweep_order.empty())
  {
    chi_log.Log(LOG_0ERROR) << "Empty linear sweep ordering.";
    exit(EXIT_FAILURE);
  }


  //================================= Compute reorder mapping
  // This mapping allows us to punch in
  // the location id and find what its
  // id is in the TDG
  size_t num_ord = glob_linear_sweep_order.size();
  for (int k=0; k<num_ord; k++)
  {
    int loc = glob_linear_sweep_order[k];
    glob_order_mapping[loc] = k;
  }

  //================================= Determine sweep order ranks
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Determining sweep order ranks.";
  int abs_max_rank = 0;
  for (int k=0; k<num_ord; k++)
  {
    chi_log.Log(LOG_0VERBOSE_1) << k;
    int loc = glob_linear_sweep_order[k];
    if (global_dependencies[loc].empty())
    {
      glob_sweep_order_rank[k] = 0;
    }
    else
    {
      int max_rank = -1;
      for (int dep=0; dep<global_dependencies[loc].size(); dep++)
      {
        int dep_loc = global_dependencies[loc][dep];
        if (dep_loc <0)
          continue;
        int dep_mapped_index = glob_order_mapping[dep_loc];

        if (glob_sweep_order_rank[dep_mapped_index] > max_rank)
        {
          max_rank = glob_sweep_order_rank[dep_mapped_index];
        }
      }
      glob_sweep_order_rank[k] = max_rank + 1;
      if ((max_rank + 1) > abs_max_rank)
      {
        abs_max_rank = max_rank + 1;
      }
    }
  }

  //================================= Generate TDG structure
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Generating TDG structure.";
  for (int r=0; r<=abs_max_rank; r++)
  {
    auto new_stdg = new chi_mesh::sweep_management::STDG;
    sweep_order->global_sweep_planes.push_back(new_stdg);

    for (int k=0; k<num_ord; k++)
    {
      if (glob_sweep_order_rank[k] == r)
      {
        new_stdg->item_id.push_back(glob_linear_sweep_order[k]);
      }
    }
  }

  TDG.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Done computing sweep ordering.";

  return sweep_order;
}
