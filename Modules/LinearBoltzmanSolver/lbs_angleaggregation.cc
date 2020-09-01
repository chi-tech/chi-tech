#include "lbs_linear_boltzman_solver.h"

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

//###################################################################
/**Initializes angle aggregation for a groupset.*/
void LinearBoltzman::Solver::InitAngleAggPolar(LBSGroupset *groupset)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Initializing angle aggregation: Polar";

  if (groupset->quadrature->type !=
    chi_math::AngularQuadratureType::ProductQuadrature)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to initialize angle aggregation. "
         "Polar angle aggregation cannot be used by the current "
         "angular quadrature.";
    exit(EXIT_FAILURE);
  }

  auto product_quadrature =
    std::static_pointer_cast<chi_math::ProductQuadrature>(groupset->quadrature);


  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int d_azi   = std::max((int)(product_quadrature->azimu_ang.size()/4),1);
  int num_azi = product_quadrature->azimu_ang.size();
  int num_pol = product_quadrature->polar_ang.size();
  int pa = num_pol/2;

  int num_angset_grps = 4; //Default Extruded and 2D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
    num_angset_grps = 1;

  //=========================================== Passing the sweep boundaries
  //                                            to the angle aggregation
  groupset->angle_agg->sim_boundaries          = sweep_boundaries;
  groupset->angle_agg->number_of_groups        = groupset->groups.size();
  groupset->angle_agg->number_of_group_subsets = groupset->grp_subsets.size();
  groupset->angle_agg->quadrature              = groupset->quadrature;
  groupset->angle_agg->grid                    = grid;

  //=========================================== Set angle aggregation
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each top hemisphere quadrant
  {
    auto angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      bool make_primary = true;
      chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

      for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
      {
        for (int an_ss=0; an_ss<groupset->ang_subsets_top.size(); an_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          AngSubSet angle_subset = groupset->ang_subsets_top[an_ss];
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
                  PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

            primary_fluds->InitializeAlphaElements(sweep_orderings[a]);
            primary_fluds->InitializeBetaElements(sweep_orderings[a]);

            fluds = primary_fluds;
          } else
          {
            fluds = new chi_mesh::sweep_management::
              AUX_FLUDS(*primary_fluds,groupset->grp_subset_sizes[gs_ss]);
          }

          auto angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a],
                          fluds,
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group->angle_sets.push_back(angleSet);
        }//for an_ss
      }//for gs_ss
    } //azi
  }//for q top
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each bot hemisphere quadrant
  {
    auto angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      bool make_primary = true;
      chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

      for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
      {
        for (int an_ss=0; an_ss<groupset->ang_subsets_bot.size(); an_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          AngSubSet angle_subset = groupset->ang_subsets_bot[an_ss];
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
            PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

            primary_fluds->InitializeAlphaElements(sweep_orderings[a+num_azi]);
            primary_fluds->InitializeBetaElements(sweep_orderings[a+num_azi]);

            fluds = primary_fluds;
          } else
          {
            fluds = new chi_mesh::sweep_management::
            AUX_FLUDS(*primary_fluds,groupset->grp_subset_sizes[gs_ss]);
          }

          auto angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a+num_azi],
                          fluds,
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group->angle_sets.push_back(angleSet);
        }//for an_ss
      }//for gs_ss
    } //azi
  }//for q bot
}

//###################################################################
/**Initializes angle aggregation for a groupset.*/
void LinearBoltzman::Solver::InitAngleAggSingle(LBSGroupset *groupset)
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Initializing angle aggregation: Single";

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  if (groupset->quadrature->type ==
      chi_math::AngularQuadratureType::ProductQuadrature)
  {
    auto product_quadrature =
      std::static_pointer_cast<chi_math::ProductQuadrature>(groupset->quadrature);

    int d_azi   = std::max((int)(product_quadrature->azimu_ang.size()/4),1);
    int num_azi = product_quadrature->azimu_ang.size();
    int num_pol = product_quadrature->polar_ang.size();
    int pa = std::max(1,num_pol/2);

    int num_angset_grps = 4; //Default Extruded and 2D
    if (options.geometry_type == GeometryType::ONED_SLAB)
        num_angset_grps = 1;

    //=========================================== Passing the sweep boundaries
    //                                            to the angle aggregation
    groupset->angle_agg->sim_boundaries          = sweep_boundaries;
    groupset->angle_agg->number_of_groups        = groupset->groups.size();
    groupset->angle_agg->number_of_group_subsets = groupset->grp_subsets.size();
    groupset->angle_agg->quadrature              = groupset->quadrature;
    groupset->angle_agg->grid                    = grid;

    //=========================================== Set angle aggregation
    for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each top hemisphere quadrant
    {
      auto angle_set_group = new TAngleSetGroup;
      groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

      for (int azi=0; azi<num_azi/num_angset_grps; azi++)
      {
        bool make_primary = true;
        chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

        for (int pr=0; pr<pa; pr++)
        {
          for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
          {
            std::vector<int> angle_indices;

            //============================================= Each quadrant gets 1/4 of azi angles
            int a = azi+q*d_azi;   //azimuthal angle index
            int p = pa - pr - 1;
            int angle_num = product_quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);

            chi_mesh::sweep_management::FLUDS* fluds;
            if (make_primary)
            {
              primary_fluds = new chi_mesh::sweep_management::
              PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

              primary_fluds->InitializeAlphaElements(sweep_orderings[angle_num]);
              primary_fluds->InitializeBetaElements(sweep_orderings[angle_num]);

              fluds = primary_fluds;
            }

            auto angleSet =
              new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                            gs_ss,
                            sweep_orderings[angle_num],
                            fluds,
                            angle_indices,
                            sweep_boundaries,
                            options.sweep_eager_limit,
                            &grid->GetCommunicator());

            angle_set_group->angle_sets.push_back(angleSet);
          }
        }//for pr

      } //azi
    }//for q top

    for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each bot hemisphere quadrant
    {
      auto angle_set_group = new TAngleSetGroup;
      groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

      for (int azi=0; azi<num_azi/num_angset_grps; azi++)
      {
        bool make_primary = true;
        chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

        for (int pr=0; pr<pa; pr++)
        {
          for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
          {
            std::vector<int> angle_indices;

            //============================================= Each quadrant gets 1/4 of azi angles
            int a = azi+q*d_azi;   //azimuthal angle index
            int p = pa + pr;
            int angle_num = product_quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);

            chi_mesh::sweep_management::FLUDS* fluds;
            if (make_primary)
            {
              primary_fluds = new chi_mesh::sweep_management::
              PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

              primary_fluds->InitializeAlphaElements(sweep_orderings[angle_num]);
              primary_fluds->InitializeBetaElements(sweep_orderings[angle_num]);

              fluds = primary_fluds;
            }

            auto angleSet =
              new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                            gs_ss,
                            sweep_orderings[angle_num],
                            fluds,
                            angle_indices,
                            sweep_boundaries,
                            options.sweep_eager_limit,
                            &grid->GetCommunicator());

            angle_set_group->angle_sets.push_back(angleSet);
          }
        }//for pr

      } //azi
    }//for q bot
  }//Product Quadrature
  else if (groupset->quadrature->type !=
           chi_math::AngularQuadratureType::ProductQuadrature)
  {
    //=========================================== Passing the sweep boundaries
    //                                            to the angle aggregation
    groupset->angle_agg->sim_boundaries          = sweep_boundaries;
    groupset->angle_agg->number_of_groups        = groupset->groups.size();
    groupset->angle_agg->number_of_group_subsets = groupset->grp_subsets.size();
    groupset->angle_agg->quadrature              = groupset->quadrature;
    groupset->angle_agg->grid                    = grid;

    //=========================================== Set angle aggregation
    for (int q=0; q<1; q++)  //%%%%%%%%% Just a single group
    {
      auto angle_set_group = new TAngleSetGroup;
      groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

      for (int n=0; n<groupset->quadrature->abscissae.size(); ++n)
      {
        bool make_primary = true;
        chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

        for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
        {
          std::vector<int> angle_indices;

          angle_indices.push_back(n);

          chi_mesh::sweep_management::FLUDS* fluds;
          if (make_primary)
          {
            primary_fluds = new chi_mesh::sweep_management::
            PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

            primary_fluds->InitializeAlphaElements(sweep_orderings[n]);
            primary_fluds->InitializeBetaElements(sweep_orderings[n]);

            fluds = primary_fluds;
          }

          auto angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[n],
                          fluds,
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group->angle_sets.push_back(angleSet);
        }

      } //angle
    }//for q top
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Failed to initialize angle aggregation. "
         "Unsupported angular quadrature.";
    exit(EXIT_FAILURE);
  }





}