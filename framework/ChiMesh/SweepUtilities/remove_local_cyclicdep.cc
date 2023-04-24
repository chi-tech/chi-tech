#include "../chi_mesh.h"

#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiGraph/chi_directed_graph.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"


//###################################################################
/**Removes local cyclic dependencies.*/
void chi_mesh::sweep_management::
  RemoveLocalCyclicDependencies(std::shared_ptr<SPDS> sweep_order,
                                chi_graph::DirectedGraph &local_DG)
{
  auto edges_to_remove = local_DG.RemoveCyclicDependencies();

  for (auto& edge_to_remove : edges_to_remove)
  {
    sweep_order->local_cyclic_dependencies.emplace_back(
      edge_to_remove.first,
      edge_to_remove.second);
  }
}