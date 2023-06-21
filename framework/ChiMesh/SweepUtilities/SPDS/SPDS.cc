#include "SPDS.h"

#include "ChiGraph/chi_directed_graph.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiConsole/chi_console.h"
#include "ChiTimer/chi_timer.h"

#include <algorithm>


//###################################################################
/** Given a location J index, maps to a predecessor location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToPrelocI(int locJ) const
{
  for (int i=0; i<location_dependencies.size(); i++)
  {
    if (location_dependencies[i] == locJ)
    {
      return i;
    }
  }

  for (int i=0; i<delayed_location_dependencies.size(); i++)
  {
    if (delayed_location_dependencies[i] == locJ)
    {
      return -(i+1);
    }
  }

  Chi::log.LogAllError()
    << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}

//###################################################################
/** Given a location J index, maps to a dependent location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToDeplocI(int locJ) const
{
  for (int i=0; i<location_successors.size(); i++)
  {
    if (location_successors[i] == locJ)
    {
      return i;
    }
  }

  Chi::log.LogAllError()
    << "SPDS Invalid mapping encountered in MapLocJToDeplocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}


//###################################################################
/**Builds the task dependency graph.*/
void chi_mesh::sweep_management::SPDS::BuildTaskDependencyGraph(bool cycle_allowance_flag)
{

  std::vector<std::pair<int,int>> edges_to_remove;
  std::vector<int> raw_edges_to_remove;
  chi_graph::DirectedGraph TDG;

  //============================================= Build graph on home location
  if (Chi::mpi.location_id == 0)
  {
    Chi::log.Log0Verbose1()
      << Chi::program_timer.GetTimeString()
      << " Building Task Dependency Graphs.";

    //====================================== Add vertices to the graph
    for (int loc=0; loc< Chi::mpi.process_count; loc++)
      TDG.AddVertex();

    //====================================== Add dependencies
    for (int loc=0; loc< Chi::mpi.process_count; loc++)
      for (int dep=0; dep<global_dependencies[loc].size(); dep++)
        TDG.AddEdge(global_dependencies[loc][dep], loc);

    //====================================== Remove cyclic dependencies
    if (cycle_allowance_flag)
    {
      Chi::log.Log0Verbose1()
        << Chi::program_timer.GetTimeString()
        << " Removing intra-cellset cycles.";
      auto edges_to_remove_temp = TDG.RemoveCyclicDependencies();
      for (const auto& [v0,v1] : edges_to_remove_temp)
        edges_to_remove.emplace_back(v0,v1);

    }

    //====================================== Serialize edges to be removed
    raw_edges_to_remove.resize(edges_to_remove.size()*2,0);
    int i=0;
    for (const auto& edge : edges_to_remove)
    {
      raw_edges_to_remove[i++] = edge.first;
      raw_edges_to_remove[i++] = edge.second;
    }
  }//if home

  //============================================= Broadcast edge buffer size
  int edge_buffer_size = 0;

  if (Chi::mpi.location_id == 0)
    edge_buffer_size = static_cast<int>(raw_edges_to_remove.size());

  MPI_Bcast(&edge_buffer_size,      //Buffer
            1, MPI_INT,             //Count and datatype
            0,                      //Root location
            MPI_COMM_WORLD);        //Communicator

  //============================================= Broadcast edges
  if (Chi::mpi.location_id != 0)
    raw_edges_to_remove.resize(edge_buffer_size,-1);

  MPI_Bcast(raw_edges_to_remove.data(),      //Buffer
            edge_buffer_size, MPI_INT, //Count and datatype
            0,                         //Root location
            MPI_COMM_WORLD);           //Communicator

  //============================================= De-serialize edges
  if (Chi::mpi.location_id != 0)
  {
    edges_to_remove.resize(edge_buffer_size/2,std::pair<int,int>(0,0));
    int i = 0;
    for (auto& edge : edges_to_remove)
    {
      edge.first  = raw_edges_to_remove[i++];
      edge.second = raw_edges_to_remove[i++];
    }
  }

  //============================================= Remove edges
  for (auto& edge_to_remove : edges_to_remove)
  {
    int rlocI  = edge_to_remove.first;
    int locI = edge_to_remove.second;

    if (Chi::mpi.location_id == 0)
      TDG.RemoveEdge(rlocI, locI);

    if (locI == Chi::mpi.location_id)
    {
      auto dependent_location =
        std::find(location_dependencies.begin(),
                  location_dependencies.end(),
                  rlocI);
      location_dependencies.erase(dependent_location);
      delayed_location_dependencies.push_back(rlocI);
    }

    if (rlocI == Chi::mpi.location_id)
      delayed_location_successors.push_back(locI);
  }

  //============================================= Generate topological sort
  std::vector<int> glob_linear_sweep_order;
  if (Chi::mpi.location_id == 0)
  {
    Chi::log.LogAllVerbose2()
      << Chi::program_timer.GetTimeString()
      << "   - Generating topological sort.";
    auto so_temp = TDG.GenerateTopologicalSort();
    for (auto v : so_temp)
      glob_linear_sweep_order.emplace_back(v);

    if (glob_linear_sweep_order.empty())
    {
      Chi::log.LogAllError()
        << "Topological sorting for global sweep-ordering failed. "
        << "Cyclic dependencies detected. Cycles need to be allowed"
        << " by calling application.";
      Chi::Exit(EXIT_FAILURE);
    }
  }

  //============================================= Broadcasting topsort size
  int topsort_buffer_size = 0;

  if (Chi::mpi.location_id == 0)
    topsort_buffer_size = glob_linear_sweep_order.size();

  MPI_Bcast(&topsort_buffer_size,   //Buffer
            1, MPI_INT,             //Count and datatype
            0,                      //Root location
            MPI_COMM_WORLD);        //Communicator

  //============================================= Broadcast topological sort
  if (Chi::mpi.location_id != 0)
    glob_linear_sweep_order.resize(topsort_buffer_size,-1);

  MPI_Bcast(glob_linear_sweep_order.data(),//Buffer
            topsort_buffer_size, MPI_INT,  //Count and datatype
            0,                             //Root location
            MPI_COMM_WORLD);               //Communicator

  //============================================= Compute reorder mapping
  // This mapping allows us to punch in
  // the location id and find what its
  // id is in the TDG
  std::vector<int> glob_order_mapping(Chi::mpi.process_count,-1);

  for (int k=0; k< Chi::mpi.process_count; k++)
  {
    int loc = glob_linear_sweep_order[k];
    glob_order_mapping[loc] = k;
  }

  //============================================= Determine sweep order ranks
  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Determining sweep order ranks.";

  std::vector<int> glob_sweep_order_rank(Chi::mpi.process_count,-1);

  int abs_max_rank = 0;
  for (int k=0; k< Chi::mpi.process_count; k++)
  {
    int loc = glob_linear_sweep_order[k];
    if (global_dependencies[loc].empty())
      glob_sweep_order_rank[k] = 0;
    else
    {
      int max_rank = -1;
      for (auto dep_loc : global_dependencies[loc])
      {
        if (dep_loc <0) continue;
        int dep_mapped_index = glob_order_mapping[dep_loc];

        if (glob_sweep_order_rank[dep_mapped_index] > max_rank)
          max_rank = glob_sweep_order_rank[dep_mapped_index];
      }
      glob_sweep_order_rank[k] = max_rank + 1;
      if ((max_rank + 1) > abs_max_rank)
        abs_max_rank = max_rank + 1;
    }
  }

  //============================================= Generate TDG structure
  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Generating TDG structure.";
  for (int r=0; r<=abs_max_rank; r++)
  {
    chi_mesh::sweep_management::STDG new_stdg;

    for (int k=0; k< Chi::mpi.process_count; k++)
    {
      if (glob_sweep_order_rank[k] == r)
        new_stdg.item_id.push_back(glob_linear_sweep_order[k]);
    }
    global_sweep_planes.push_back(new_stdg);
  }
}