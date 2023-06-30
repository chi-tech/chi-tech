#include "lbs_discrete_ordinates_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Tools/lbs_make_subset.h"

#include "mesh/SweepUtilities/FLUDS/AUX_FLUDS.h"

#include "console/chi_console.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <iomanip>

#include "utils/chi_timer.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

//###################################################################
/**Initializes fluds data structures.*/

void lbs::DiscreteOrdinatesSolver::
  InitFluxDataStructures(LBSGroupset& groupset)
{
  namespace sweep_namespace = chi_mesh::sweep_management;
  typedef sweep_namespace::AngleSetGroup TAngleSetGroup;
  typedef sweep_namespace::AngleSet TAngleSet;

  const auto& quadrature_sweep_info =
    quadrature_unq_so_grouping_map_[groupset.quadrature_];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map    = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.groups_.size();
  const size_t gs_num_ss = groupset.grp_subset_infos_.size();

  //=========================================== Passing the sweep boundaries
  //                                            to the angle aggregation
  typedef chi_mesh::sweep_management::AngleAggregation AngleAgg;
  groupset.angle_agg_ = std::make_shared<AngleAgg>(sweep_boundaries_,
                                                   gs_num_grps,
                                                   gs_num_ss,
                                                   groupset.quadrature_,
                                                   grid_ptr_);

  TAngleSetGroup angle_set_group;
  for (const auto& so_grouping : unique_so_groupings)
  {
    const size_t master_dir_id = so_grouping.front();
    const size_t so_id = dir_id_to_so_map.at(master_dir_id);

    const auto& sweep_ordering =
      quadrature_spds_map_[groupset.quadrature_][so_id];
    const auto& fluds_template =
      *quadrature_fluds_templates_map_[groupset.quadrature_][so_id];

    //Compute direction subsets
    const auto dir_subsets = lbs::MakeSubSets(so_grouping.size(),
                                              groupset.master_num_ang_subsets_);

    for (size_t gs_ss=0; gs_ss<gs_num_ss; gs_ss++)
    {
      const size_t gs_ss_size = groupset.grp_subset_infos_[gs_ss].ss_size;
      for (const auto & dir_ss_info : dir_subsets)
      {
        const auto& dir_ss_begin = dir_ss_info.ss_begin;
        const auto& dir_ss_end   = dir_ss_info.ss_end;
        const auto& dir_ss_size   = dir_ss_info.ss_size;

        std::vector<size_t> angle_indices(dir_ss_size, 0);
        {
          size_t k = 0;
          for (size_t n=dir_ss_begin; n<=dir_ss_end; ++n)
            angle_indices[k++] = so_grouping[n];
        }

        using namespace chi_mesh::sweep_management;
        auto aux_fluds = std::make_shared<AUX_FLUDS>(fluds_template,gs_ss_size);

        auto fluds = std::dynamic_pointer_cast<FLUDS>(aux_fluds);

        if (not fluds)
          throw std::runtime_error(std::string(__PRETTY_FUNCTION__) +
          ": Casting failure.");

        auto angleSet = std::make_shared<TAngleSet>(
          gs_ss_size,
          gs_ss,
          *sweep_ordering,
          fluds,
          angle_indices,
          sweep_boundaries_,
          options_.sweep_eager_limit,
          *grid_local_comm_set_);

        angle_set_group.angle_sets.push_back(angleSet);
      }//for an_ss
    }//for gs_ss
  }//for so_grouping

  groupset.angle_agg_->angle_set_groups.push_back(std::move(angle_set_group));

  if (options_.verbose_inner_iterations)
    Chi::log.Log()
      << Chi::program_timer.GetTimeString()
      << " Initialized Angle Aggregation.   "
      << "         Process memory = "
      << std::setprecision(3) << chi::ChiConsole::GetMemoryUsageInMB()
      << " MB.";


  Chi::mpi.Barrier();
}
