#include "../chi_mesh.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiGraph/chi_directed_graph.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <algorithm>

//###################################################################
/**Removes local cyclic dependencies.*/
void chi_mesh::sweep_management::
  RemoveLocalCyclicDependencies(chi_mesh::sweep_management::SPDS *sweep_order,
                                chi_graph::DirectedGraph &local_DG)
{
  //============================================= Find initial SCCs
  auto SCCs = local_DG.FindStronglyConnectedConnectionns();

  //============================================= Rinse remove edges
  std::vector<std::pair<int,int>> edges_to_remove;
  int iter=0;
  while (not SCCs.empty())
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "Inter cell cyclic dependency removal. Iteration " << ++iter;
    //=================================== Loop over sub-graphs
    edges_to_remove.clear();
    for (auto& subDG : SCCs)
    {
      //Identify edges to remove
      bool found_edge_to_remove = false;
      for (int v : subDG)
      {
        for (int w : local_DG.vertices[v].ds_edge)
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
      local_DG.RemoveEdge(edge_to_remove.first, edge_to_remove.second);
      sweep_order->local_cyclic_dependencies.emplace_back(
        edge_to_remove.first,
        edge_to_remove.second);
    }

    // Refind SCCs
    SCCs = local_DG.FindStronglyConnectedConnectionns();
  }
}

//###################################################################
/**Removes global cyclic dependencies.*/
void chi_mesh::sweep_management::
  RemoveGlobalCyclicDependencies(
    chi_mesh::sweep_management::SPDS* sweep_order,
    chi_graph::DirectedGraph& TDG)
{
  //============================================= Find initial SCCs
  auto SCCs = TDG.FindStronglyConnectedConnectionns();

  //============================================= Rinse remove edges
  std::vector<std::pair<int,int>> edges_to_remove;
  int iter=0;
  while (not SCCs.empty())
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "Inter cell-set cyclic dependency removal. Iteration " << ++iter;
    //=================================== Loop over sub-graphs
    edges_to_remove.clear();
    for (auto& subDG : SCCs)
    {
      //Identify edges to remove
      bool found_edge_to_remove = false;
      for (int v : subDG)
      {
        for (int w : TDG.vertices[v].ds_edge)
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
      int rlocI  = edge_to_remove.first;
      int locI = edge_to_remove.second;
      TDG.RemoveEdge(rlocI, locI);

      if (locI == chi_mpi.location_id)
      {
        auto dependent_location =
          std::find(sweep_order->location_dependencies.begin(),
                    sweep_order->location_dependencies.end(),
                    rlocI);
        sweep_order->location_dependencies.erase(dependent_location);
        sweep_order->delayed_location_dependencies.push_back(rlocI);
      }

      if (rlocI == chi_mpi.location_id)
      {
        sweep_order->delayed_location_successors.push_back(locI);
      }
    }

    // Refind SCCs
    SCCs = TDG.FindStronglyConnectedConnectionns();
  }
}