#include "lbs_linear_boltzman_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>

#include <algorithm>

#include <chi_log.h>

extern ChiLog chi_log;

typedef chi_mesh::SweepManagement::AngleSet TAngleSet;
typedef chi_mesh::SweepManagement::AngleSetGroup TAngleSetGroup;



//###################################################################
/**Initializes angle aggregation for a groupset.*/
void LinearBoltzman::Solver::InitAngleAggPolar(LBSGroupset *groupset)
{
  chi_log.Log(LOG_0) << "Initializing angle aggregation: Polar";

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int d_azi   = std::max((int)(groupset->quadrature->azimu_ang.size()/4),1);
  int num_azi = groupset->quadrature->azimu_ang.size();
  int num_pol = groupset->quadrature->polar_ang.size();
  int pa = num_pol/2;

  int num_angset_grps = 4; //Default Extruded and 2D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
    num_angset_grps = 1;




  //=========================================== Set angle aggregation
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each top hemisphere quadrant
  {
    TAngleSetGroup* angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
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
            int angle_num = groupset->quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);
          }//for pr


          TAngleSet* angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a],
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &comm_set);

          angle_set_group->angle_sets.push_back(angleSet);
        }//for an_ss
      }//for gs_ss
    } //azi
  }//for q top
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each bot hemisphere quadrant
  {
    TAngleSetGroup* angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
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
            int angle_num = groupset->quadrature->GetAngleNum(p,a);
            angle_indices.push_back(angle_num);
          }//for pr
          TAngleSet* angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a+num_azi],
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &comm_set);

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
  chi_log.Log(LOG_0) << "Initializing angle aggregation: Single";

  int d_azi   = std::max((int)(groupset->quadrature->azimu_ang.size()/4),1);
  int num_azi = groupset->quadrature->azimu_ang.size();
  int num_pol = groupset->quadrature->polar_ang.size();
  int pa = std::max(1,num_pol/2);

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int num_angset_grps = 4; //Default Extruded and 2D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
    num_angset_grps = 1;

  //=========================================== Set angle aggregation
  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each top hemisphere quadrant
  {
    TAngleSetGroup* angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      for (int pr=0; pr<pa; pr++)
      {
        for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          int p = pa - pr - 1;
          int angle_num = groupset->quadrature->GetAngleNum(p,a);
          angle_indices.push_back(angle_num);
          TAngleSet* angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a],
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &comm_set);

          angle_set_group->angle_sets.push_back(angleSet);
        }
      }//for pr

    } //azi
  }//for q top

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
    return;

  for (int q=0; q<num_angset_grps; q++)  //%%%%%%%%% for each bot hemisphere quadrant
  {
    TAngleSetGroup* angle_set_group = new TAngleSetGroup;
    groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

    for (int azi=0; azi<num_azi/num_angset_grps; azi++)
    {
      for (int pr=0; pr<pa; pr++)
      {
        for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
        {
          std::vector<int> angle_indices;

          //============================================= Each quadrant gets 1/4 of azi angles
          int a = azi+q*d_azi;   //azimuthal angle index
          int p = pa + pr;
          int angle_num = groupset->quadrature->GetAngleNum(p,a);
          angle_indices.push_back(angle_num);
          TAngleSet* angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings[a+num_azi],
                          angle_indices,
                          sweep_boundaries,
                          options.sweep_eager_limit,
                          &comm_set);

          angle_set_group->angle_sets.push_back(angleSet);
        }
      }//for pr

    } //azi
  }//for q bot
}