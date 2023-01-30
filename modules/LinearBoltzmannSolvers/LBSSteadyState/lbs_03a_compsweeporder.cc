#include "lbs_linear_boltzmann_solver.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

typedef chi_mesh::sweep_management::AngleSet TAngleSet;
typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;

#include "ChiConsole/chi_console.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"


#include <iomanip>

//###################################################################
/**Initializes the sweep ordering for the given groupset.*/
void lbs::SteadyStateSolver::ComputeSweepOrderings(LBSGroupset& groupset) const
{
  if (options.verbose_inner_iterations)
    chi::log.Log()
      << chi::program_timer.GetTimeString()
      << " Computing Sweep ordering.\n";

  const auto unq_groupings_and_mapping =
    AssociateSOsAndDirections(*grid,
                              *groupset.quadrature,
                              groupset.angleagg_method,
                              options.geometry_type);
  const auto& unique_so_groupings = unq_groupings_and_mapping.first;
  const auto& dir_id_to_so_map = unq_groupings_and_mapping.second;

  groupset.unique_so_groupings = unique_so_groupings;
  groupset.dir_id_to_so_map    = dir_id_to_so_map;

  //============================================= Clear sweep ordering
  groupset.sweep_orderings.clear();
  groupset.sweep_orderings.shrink_to_fit();

  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  auto mesher = mesh_handler.volume_mesher;

  bool no_cycles_parmetis_partitioning =
    (mesher->options.partition_type ==
     chi_mesh::VolumeMesher::PartitionType::PARMETIS
     and (not groupset.allow_cycles));

  bool is_1D_geometry = options.geometry_type == GeometryType::ONED_SLAB;

  //============================================= Check possibility of cycles
  if (no_cycles_parmetis_partitioning and
      not is_1D_geometry and
      chi::mpi.process_count>1)
  {
    chi::log.LogAllError()
      << "When using PARMETIS type partitioning then groupset iterative method"
         " must be NPT_CLASSICRICHARDSON_CYCLES or NPT_GMRES_CYCLES";
    chi::Exit(EXIT_FAILURE);
  }

  //================================================== Compute sweep orderings
  for (const auto& so_grouping : unique_so_groupings)
  {
    if (so_grouping.empty()) continue;

    const size_t master_dir_id = so_grouping.front();
    const auto& omega = groupset.quadrature->omegas[master_dir_id];
    const auto new_swp_order =
      chi_mesh::sweep_management::
      CreateSweepOrder(omega,
                       this->grid,
                       groupset.allow_cycles);
    groupset.sweep_orderings.emplace_back(new_swp_order);
  }

  if (options.verbose_inner_iterations)
    chi::log.Log()
      << chi::program_timer.GetTimeString()
      << " Done computing sweep orderings.           Process memory = "
      << std::setprecision(3)
      << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";

}