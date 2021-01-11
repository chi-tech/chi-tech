#include "SPDS.h"

#include "ChiGraph/chi_directed_graph.h"

#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiConsole/chi_console.h"
#include "ChiTimer/chi_timer.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiConsole&  chi_console;
extern ChiTimer   chi_program_timer;

#include <algorithm>


//###################################################################
/** Given a location J index, maps to a predecessor location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToPrelocI(int locJ)
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

  chi_log.Log(LOG_ALLERROR)
    << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  exit(EXIT_FAILURE);
}

//###################################################################
/** Given a location J index, maps to a dependent location.*/
int chi_mesh::sweep_management::SPDS::MapLocJToDeplocI(int locJ)
{
  for (int i=0; i<location_successors.size(); i++)
  {
    if (location_successors[i] == locJ)
    {
      return i;
    }
  }

  chi_log.Log(LOG_ALLERROR)
    << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  exit(EXIT_FAILURE);
  return -1;
}

////###################################################################
///** Adds a location dependency for this location.*/
//void chi_mesh::sweep_management::SPDS::AddLocalDependecy(int location_index)
//{
//  if (location_index<0){return;}
//
//  bool already_there = false;
//  size_t num_deps = location_dependencies.size();
//  for (int k=0; k<num_deps; k++)
//  {
//    if (location_dependencies[k] == location_index)
//    {
//      already_there = true;
//      break;
//    }
//  }
//  if (!already_there)
//  {
//    location_dependencies.push_back(location_index);
//  }
//}
//
////###################################################################
///** Adds a location successor for this location.*/
//void chi_mesh::sweep_management::SPDS::AddLocalSuccessor(int location_index)
//{
//  if (location_index<0){return;}
//
//  bool already_there = false;
//  size_t num_sucs = location_successors.size();
//  for (int k=0; k<num_sucs; k++)
//  {
//    if (location_successors[k] == location_index)
//    {
//      already_there = true;
//      break;
//    }
//  }
//  if (!already_there)
//  {
//    location_successors.push_back(location_index);
//  }
//}

//###################################################################
/**Builds the task dependency graph.*/
void chi_mesh::sweep_management::SPDS::BuildTaskDependencyGraph(bool cycle_allowance_flag)
{
//  chi_log.Log(LOG_0VERBOSE_1)
//    << chi_program_timer.GetTimeString()
//    << " Building Task Dependency Graphs.";
//  chi_graph::DirectedGraph TDG;
//
//  //============================================= Add vertices to the graph
//  for (int loc=0; loc<chi_mpi.process_count; loc++)
//    TDG.AddVertex();
//
//  //============================================= Add dependencies
//  for (int loc=0; loc<chi_mpi.process_count; loc++)
//    for (int dep=0; dep<global_dependencies[loc].size(); dep++)
//      TDG.AddEdge(global_dependencies[loc][dep], loc);
//
//  //============================================= Filter dependencies
//  //                                              for cycles
//  if (cycle_allowance_flag)
//  {
//    chi_log.Log(LOG_0VERBOSE_1)
//      << chi_program_timer.GetTimeString()
//      << " Removing intra-cellset cycles.";
//    RemoveGlobalCyclicDependencies(this,TDG);
//  }
//
//  //============================================= Generate topological sort
//  chi_log.Log(LOG_ALLVERBOSE_2)
//    << chi_program_timer.GetTimeString()
//    << "   - Generating topological sort.";
//  std::vector<int> glob_linear_sweep_order = TDG.GenerateTopologicalSort();
//
//  if (glob_linear_sweep_order.empty())
//  {
//    chi_log.Log(LOG_ALLERROR)
//      << "Topological sorting for global sweep-ordering failed. "
//      << "Cyclic dependencies detected. Cycles need to be allowed"
//      << " by calling application.";
//    exit(EXIT_FAILURE);
//  }
//
//  //============================================= Compute reorder mapping
//  // This mapping allows us to punch in
//  // the location id and find what its
//  // id is in the TDG
//  std::vector<int> glob_order_mapping(chi_mpi.process_count,-1);
//
//  for (int k=0; k<chi_mpi.process_count; k++)
//  {
//    int loc = glob_linear_sweep_order[k];
//    glob_order_mapping[loc] = k;
//  }
//
//  //============================================= Determine sweep order ranks
//  chi_log.Log(LOG_0VERBOSE_1)
//    << chi_program_timer.GetTimeString()
//    << " Determining sweep order ranks.";
//
//  std::vector<int> glob_sweep_order_rank(chi_mpi.process_count,-1);
//
//  int abs_max_rank = 0;
//  for (int k=0; k<chi_mpi.process_count; k++)
//  {
//    int loc = glob_linear_sweep_order[k];
//    if (global_dependencies[loc].empty())
//      glob_sweep_order_rank[k] = 0;
//    else
//    {
//      int max_rank = -1;
//      for (auto dep_loc : global_dependencies[loc])
//      {
//        if (dep_loc <0) continue;
//        int dep_mapped_index = glob_order_mapping[dep_loc];
//
//        if (glob_sweep_order_rank[dep_mapped_index] > max_rank)
//          max_rank = glob_sweep_order_rank[dep_mapped_index];
//      }
//      glob_sweep_order_rank[k] = max_rank + 1;
//      if ((max_rank + 1) > abs_max_rank)
//        abs_max_rank = max_rank + 1;
//    }
//  }
//
//  //============================================= Generate TDG structure
//  chi_log.Log(LOG_0VERBOSE_1)
//    << chi_program_timer.GetTimeString()
//    << " Generating TDG structure.";
//  for (int r=0; r<=abs_max_rank; r++)
//  {
//    auto new_stdg = new chi_mesh::sweep_management::STDG;
//    global_sweep_planes.push_back(new_stdg);
//
//    for (int k=0; k<chi_mpi.process_count; k++)
//    {
//      if (glob_sweep_order_rank[k] == r)
//        new_stdg->item_id.push_back(glob_linear_sweep_order[k]);
//    }
//  }

  std::vector<std::pair<int,int>> edges_to_remove;
  std::vector<int> raw_edges_to_remove;
  chi_graph::DirectedGraph TDG;

  //============================================= Build graph on home location
  if (chi_mpi.location_id == 0)
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << chi_program_timer.GetTimeString()
      << " Building Task Dependency Graphs.";

    //====================================== Add vertices to the graph
    for (int loc=0; loc<chi_mpi.process_count; loc++)
      TDG.AddVertex();

    //====================================== Add dependencies
    for (int loc=0; loc<chi_mpi.process_count; loc++)
      for (int dep=0; dep<global_dependencies[loc].size(); dep++)
        TDG.AddEdge(global_dependencies[loc][dep], loc);

    //====================================== Remove cyclic dependencies
    if (cycle_allowance_flag)
    {
      chi_log.Log(LOG_0VERBOSE_1)
        << chi_program_timer.GetTimeString()
        << " Removing intra-cellset cycles.";
      edges_to_remove = TDG.RemoveCyclicDependencies();
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

  if (chi_mpi.location_id == 0)
    edge_buffer_size = raw_edges_to_remove.size();

  MPI_Bcast(&edge_buffer_size,      //Buffer
            1, MPI_INT,             //Count and datatype
            0,                      //Root location
            MPI_COMM_WORLD);        //Communicator

  //============================================= Broadcast edges
  if (chi_mpi.location_id != 0)
    raw_edges_to_remove.resize(edge_buffer_size,-1);

  MPI_Bcast(raw_edges_to_remove.data(),      //Buffer
            edge_buffer_size, MPI_INT, //Count and datatype
            0,                         //Root location
            MPI_COMM_WORLD);           //Communicator

  //============================================= De-serialize edges
  if (chi_mpi.location_id != 0)
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

    if (chi_mpi.location_id == 0)
      TDG.RemoveEdge(rlocI, locI);

    if (locI == chi_mpi.location_id)
    {
      auto dependent_location =
        std::find(location_dependencies.begin(),
                  location_dependencies.end(),
                  rlocI);
      location_dependencies.erase(dependent_location);
      delayed_location_dependencies.push_back(rlocI);
    }

    if (rlocI == chi_mpi.location_id)
      delayed_location_successors.push_back(locI);
  }

  //============================================= Generate topological sort
  std::vector<int> glob_linear_sweep_order;
  if (chi_mpi.location_id == 0)
  {
    chi_log.Log(LOG_ALLVERBOSE_2)
      << chi_program_timer.GetTimeString()
      << "   - Generating topological sort.";
    glob_linear_sweep_order = TDG.GenerateTopologicalSort();

    if (glob_linear_sweep_order.empty())
    {
      chi_log.Log(LOG_ALLERROR)
        << "Topological sorting for global sweep-ordering failed. "
        << "Cyclic dependencies detected. Cycles need to be allowed"
        << " by calling application.";
      exit(EXIT_FAILURE);
    }
  }

  //============================================= Broadcasting topsort size
  int topsort_buffer_size = 0;

  if (chi_mpi.location_id == 0)
    topsort_buffer_size = glob_linear_sweep_order.size();

  MPI_Bcast(&topsort_buffer_size,   //Buffer
            1, MPI_INT,             //Count and datatype
            0,                      //Root location
            MPI_COMM_WORLD);        //Communicator

  //============================================= Broadcast topological sort
  if (chi_mpi.location_id != 0)
    glob_linear_sweep_order.resize(topsort_buffer_size,-1);

  MPI_Bcast(glob_linear_sweep_order.data(),//Buffer
            topsort_buffer_size, MPI_INT,  //Count and datatype
            0,                             //Root location
            MPI_COMM_WORLD);               //Communicator

  //============================================= Compute reorder mapping
  // This mapping allows us to punch in
  // the location id and find what its
  // id is in the TDG
  std::vector<int> glob_order_mapping(chi_mpi.process_count,-1);

  for (int k=0; k<chi_mpi.process_count; k++)
  {
    int loc = glob_linear_sweep_order[k];
    glob_order_mapping[loc] = k;
  }

  //============================================= Determine sweep order ranks
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Determining sweep order ranks.";

  std::vector<int> glob_sweep_order_rank(chi_mpi.process_count,-1);

  int abs_max_rank = 0;
  for (int k=0; k<chi_mpi.process_count; k++)
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
  chi_log.Log(LOG_0VERBOSE_1)
    << chi_program_timer.GetTimeString()
    << " Generating TDG structure.";
  for (int r=0; r<=abs_max_rank; r++)
  {
    chi_mesh::sweep_management::STDG new_stdg;

    for (int k=0; k<chi_mpi.process_count; k++)
    {
      if (glob_sweep_order_rank[k] == r)
        new_stdg.item_id.push_back(glob_linear_sweep_order[k]);
    }
    global_sweep_planes.push_back(new_stdg);
  }
}