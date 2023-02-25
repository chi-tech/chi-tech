#include "lbs_DO_steady_state.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#define ParallelParmetisNeedsCycles \
"When using PARMETIS type partitioning then groupset iterative method" \
" must be NPT_CLASSICRICHARDSON_CYCLES or NPT_GMRES_CYCLES"

#define IsParallel chi::mpi.process_count>1

#define IsPartitionTypeParmetis \
mesher->options.partition_type == chi_mesh::VolumeMesher::PartitionType::PARMETIS

namespace lbs
{

void DiscOrdSteadyStateSolver::InitializeSweepDataStructures()
{
  if (options_.verbose_inner_iterations)
    chi::log.Log() << chi::program_timer.GetTimeString()
                   << " Initializing sweep datastructures.\n";

  //=================================== Perform checks
  {
    auto& mesh_handler = chi_mesh::GetCurrentHandler();
    auto mesher = mesh_handler.volume_mesher;

    for (const auto& groupset : groupsets_)
    {
      bool no_cycles_parmetis_partitioning =
        (IsPartitionTypeParmetis and (not groupset.allow_cycles_));

      bool is_1D_geometry = options_.geometry_type == GeometryType::ONED_SLAB;

      if (no_cycles_parmetis_partitioning and not is_1D_geometry and IsParallel)
        throw std::logic_error(ParallelParmetisNeedsCycles);
    }//for groupset
  }

  //=================================== Define sweep ordering groups
  quadrature_unq_so_grouping_map_.clear();
  std::map<AngQuadPtr, bool> quadrature_allow_cycles_map_;
  for (auto& groupset : groupsets_)
  {
    if (quadrature_unq_so_grouping_map_.count(groupset.quadrature_) == 0)
      quadrature_unq_so_grouping_map_[groupset.quadrature_] =
        AssociateSOsAndDirections(*grid_ptr_,
                                  *groupset.quadrature_,
                                  groupset.angleagg_method_,
                                  options_.geometry_type);

    if (quadrature_allow_cycles_map_.count(groupset.quadrature_) == 0)
      quadrature_allow_cycles_map_[groupset.quadrature_] = groupset.allow_cycles_;
  }

  //=================================== Build sweep orderings
  quadrature_spds_map_.clear();
  for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
  {
    const auto& unique_so_groupings = info.first;

    for (const auto& so_grouping : unique_so_groupings)
    {
      if (so_grouping.empty()) continue;

      const size_t master_dir_id = so_grouping.front();
      const auto& omega = quadrature->omegas[master_dir_id];
      const auto new_swp_order =
        chi_mesh::sweep_management::
        CreateSweepOrder(omega,
                         this->grid_ptr_,
                         quadrature_allow_cycles_map_[quadrature]);
      quadrature_spds_map_[quadrature].push_back(new_swp_order);
    }
  }//quadrature info-pack

  //=================================== Build FLUDS templates
  quadrature_fluds_templates_map_.clear();
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
      quadrature_fluds_templates_map_[quadrature].push_back(
        std::make_shared<FLUDSTemplate>(1, grid_nodal_mappings_, *spds)
        );
  }//for quadrature spds-list pair
}

}//namespace lbs