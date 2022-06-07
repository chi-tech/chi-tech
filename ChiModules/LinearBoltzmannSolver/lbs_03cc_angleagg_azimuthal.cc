#include "lbs_linear_boltzmann_solver.h"

#include "ChiMesh/SweepUtilities/FLUDS/AUX_FLUDS.h"

#include "ChiMath/Quadratures/product_quadrature.h"

#include "chi_runtime.h"
#include "chi_log.h"

;

#include "ChiTimer/chi_timer.h"



#include "ChiConsole/chi_console.h"
#include "Groupset/lbs_groupset.h"

typedef chi_mesh::sweep_management::AngleSet TAngleSet;
typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;

#include <iomanip>

//###################################################################
/**Initializes angle aggregation for a groupset.*/
void lbs::SteadySolver::InitAngleAggAzimuthal(LBSGroupset& groupset)
{
  if (options.verbose_inner_iterations)
    chi::log.Log()
      << chi::program_timer.GetTimeString()
      << " Initializing angle aggregation: Azimuthal";

  if (groupset.quadrature->type !=
    chi_math::AngularQuadratureType::ProductQuadrature)
  {
    chi::log.LogAllError()
      << "Failed to initialize angle aggregation. "
         "Azimuthal angle aggregation cannot be used by the current "
         "angular quadrature.";
    chi::Exit(EXIT_FAILURE);
  }

  auto product_quadrature =
    std::static_pointer_cast<chi_math::ProductQuadrature>(groupset.quadrature);

  //=========================================== Passing the sweep boundaries
  //                                            to the angle aggregation
  groupset.angle_agg.Setup(sweep_boundaries,
                           groupset.groups.size(),
                           groupset.grp_subsets.size(),
                           groupset.quadrature,
                           grid);

  //=========================================== Set angle aggregation
  TAngleSetGroup angle_set_group;

  int angle_num = -1;
  for (const auto& dir_set : product_quadrature->GetDirectionMap())
  {
    const auto& n_azimu = dir_set.second.size();
    for (unsigned int quad = 0; quad < 2; ++quad)
    {
      angle_num++;

      std::vector<int> angle_indices;
      for (unsigned int n = 0; n < n_azimu/2; ++n)
        angle_indices.emplace_back(dir_set.second[quad*n_azimu/2+n]);

      bool make_primary = true;
      chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

      for (size_t gs_ss = 0; gs_ss < groupset.grp_subsets.size(); ++gs_ss)
      {
        chi_mesh::sweep_management::FLUDS* fluds;
        if (make_primary)
        {
          make_primary = false;
          primary_fluds = new chi_mesh::sweep_management::
            PRIMARY_FLUDS(groupset.grp_subset_sizes[gs_ss],
                          grid_nodal_mappings);

          chi::log.Log0Verbose1()
            << "Initializing FLUDS for omega="
            << groupset.sweep_orderings[angle_num]->omega.PrintS()
            << "         Process memory = "
            << std::setprecision(3) << chi::console.GetMemoryUsageInMB()
            << " MB.";

          primary_fluds->InitializeAlphaElements(groupset.sweep_orderings[angle_num]);
          primary_fluds->InitializeBetaElements(groupset.sweep_orderings[angle_num]);

          fluds = primary_fluds;
        }
        else
        {
          fluds = new chi_mesh::sweep_management::
          AUX_FLUDS(*primary_fluds,groupset.grp_subset_sizes[gs_ss]);
        }

        auto angleSet = std::make_shared<TAngleSet>(
                        groupset.grp_subset_sizes[gs_ss],
                        gs_ss,
                        groupset.sweep_orderings[angle_num],
                        fluds,
                        angle_indices,
                        sweep_boundaries,
                        options.sweep_eager_limit,
                        &grid->GetCommunicator());

        angle_set_group.angle_sets.push_back(angleSet);
      }
    }
  }

  groupset.angle_agg.angle_set_groups.push_back(angle_set_group);
}
