#include "../chi_mesh.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiGraph/chi_directed_graph.h"

#include "chi_runtime.h"

#include <algorithm>

//###################################################################
/**Removes global cyclic dependencies.*/
void chi_mesh::sweep_management::
  RemoveGlobalCyclicDependencies(
    chi_mesh::sweep_management::SPDS* sweep_order,
    chi_graph::DirectedGraph& TDG)
{
  auto edges_to_remove = TDG.RemoveCyclicDependencies();

  //Remove the edges
  for (auto& edge_to_remove : edges_to_remove)
  {
    int rlocI  = edge_to_remove.first;
    int locI = edge_to_remove.second;
    TDG.RemoveEdge(rlocI, locI);

    if (locI == Chi::mpi.location_id)
    {
      auto dependent_location =
        std::find(sweep_order->location_dependencies.begin(),
                  sweep_order->location_dependencies.end(),
                  rlocI);
      sweep_order->location_dependencies.erase(dependent_location);
      sweep_order->delayed_location_dependencies.push_back(rlocI);
    }

    if (rlocI == Chi::mpi.location_id)
    {
      sweep_order->delayed_location_successors.push_back(locI);
    }
  }


}