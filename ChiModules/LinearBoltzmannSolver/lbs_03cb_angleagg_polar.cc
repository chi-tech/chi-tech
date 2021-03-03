#include "lbs_linear_boltzmann_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/SweepUtilities/FLUDS/AUX_FLUDS.h>

#include "ChiMath/Quadratures/product_quadrature.h"

#include <chi_log.h>

extern ChiLog& chi_log;

typedef chi_mesh::sweep_management::AngleSet TAngleSet;
typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;

#include "ChiTimer/chi_timer.h"

extern ChiTimer chi_program_timer;

#include "ChiConsole/chi_console.h"

extern ChiConsole& chi_console;

#include "iomanip"

//###################################################################
/**Initializes angle aggregation for a groupset.*/
void LinearBoltzmann::Solver::InitAngleAggPolar(LBSGroupset& groupset)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Initializing angle aggregation: Polar";

  if (groupset.quadrature->type !=
    chi_math::AngularQuadratureType::ProductQuadrature)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to initialize angle aggregation. "
         "Polar angle aggregation cannot be used by the current "
         "angular quadrature.";
    exit(EXIT_FAILURE);
  }

  auto product_quadrature =
    std::static_pointer_cast<chi_math::ProductQuadrature>(groupset.quadrature);


  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int d_azi   = std::max((int)(product_quadrature->azimu_ang.size()/4),1);
  int num_azi = product_quadrature->azimu_ang.size();
  int num_pol = product_quadrature->polar_ang.size();
  int pa = num_pol/2;

  int num_angset_grps = 4; //Default Extruded and 2D
  if (options.geometry_type == GeometryType::ONED_SLAB)
    num_angset_grps = 1;

  //=========================================== Passing the sweep boundaries
  //                                            to the angle aggregation
  groupset.angle_agg.Setup(sweep_boundaries,
                            groupset.groups.size(),
                            groupset.grp_subsets.size(),
                            groupset.quadrature,
                            grid);

  //=========================================== Set angle aggregation
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each top hemisphere quadrant
  {
    TAngleSetGroup angle_set_group;

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      bool make_primary = true;
      chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

      for (int gs_ss=0; gs_ss<groupset.grp_subsets.size(); gs_ss++)
      {
        for (int an_ss=0; an_ss<groupset.ang_subsets_top.size(); an_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          AngSubSet angle_subset = groupset.ang_subsets_top[an_ss];
          for (int p=angle_subset.first; p<=angle_subset.second; p++)
          {
            int angle_num = product_quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);
          }//for pr

          chi_mesh::sweep_management::FLUDS* fluds;
          if (make_primary)
          {
            make_primary = false;
            primary_fluds = new chi_mesh::sweep_management::
                  PRIMARY_FLUDS(groupset.grp_subset_sizes[gs_ss],
                                grid_nodal_mappings);

            chi_log.Log(LOG_0VERBOSE_1)
              << "Initializing FLUDS for omega="
              << groupset.sweep_orderings[a]->omega.PrintS()
              << "         Process memory = "
              << std::setprecision(3) << chi_console.GetMemoryUsageInMB()
              << " MB.";

            primary_fluds->InitializeAlphaElements(groupset.sweep_orderings[a]);
            primary_fluds->InitializeBetaElements(groupset.sweep_orderings[a]);

            fluds = primary_fluds;
          } else
          {
            fluds = new chi_mesh::sweep_management::
              AUX_FLUDS(*primary_fluds,groupset.grp_subset_sizes[gs_ss]);
          }

          auto angleSet = std::make_shared<TAngleSet>(
                          groupset.grp_subset_sizes[gs_ss],
                          gs_ss,
                          groupset.sweep_orderings[a],
                          fluds,
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group.angle_sets.push_back(angleSet);
        }//for an_ss
      }//for gs_ss
    } //azi

    groupset.angle_agg.angle_set_groups.push_back(angle_set_group);
  }//for q top
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each bot hemisphere quadrant
  {
    TAngleSetGroup angle_set_group;

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      bool make_primary = true;
      chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

      for (int gs_ss=0; gs_ss<groupset.grp_subsets.size(); gs_ss++)
      {
        for (int an_ss=0; an_ss<groupset.ang_subsets_bot.size(); an_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          AngSubSet angle_subset = groupset.ang_subsets_bot[an_ss];
          for (int p=angle_subset.first; p<=angle_subset.second; p++)
          {
            int angle_num = product_quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);
          }//for pr

          chi_mesh::sweep_management::FLUDS* fluds;
          if (make_primary)
          {
            make_primary = false;
            primary_fluds = new chi_mesh::sweep_management::
            PRIMARY_FLUDS(groupset.grp_subset_sizes[gs_ss],
                          grid_nodal_mappings);

            chi_log.Log(LOG_0VERBOSE_1)
              << "Initializing FLUDS for omega="
              << groupset.sweep_orderings[a+num_azi]->omega.PrintS()
              << "         Process memory = "
              << std::setprecision(3) << chi_console.GetMemoryUsageInMB()
              << " MB.";

            primary_fluds->InitializeAlphaElements(groupset.sweep_orderings[a+num_azi]);
            primary_fluds->InitializeBetaElements(groupset.sweep_orderings[a+num_azi]);

            fluds = primary_fluds;
          } else
          {
            fluds = new chi_mesh::sweep_management::
            AUX_FLUDS(*primary_fluds,groupset.grp_subset_sizes[gs_ss]);
          }

          auto angleSet = std::make_shared<TAngleSet>(
                          groupset.grp_subset_sizes[gs_ss],
                          gs_ss,
                          groupset.sweep_orderings[a+num_azi],
                          fluds,
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group.angle_sets.push_back(angleSet);
        }//for an_ss
      }//for gs_ss
    } //azi

    groupset.angle_agg.angle_set_groups.push_back(angle_set_group);
  }//for q bot
}
