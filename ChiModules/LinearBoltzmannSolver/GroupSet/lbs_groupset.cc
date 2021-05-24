#include "lbs_groupset.h"

#include "ChiMath/Quadratures/product_quadrature.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <fstream>

//##############################################
/**Groupset constructor.*/
LBSGroupset::LBSGroupset()
{
  quadrature = nullptr;
  iterative_method = LinearBoltzmann::IterativeMethod::NONE;
  angleagg_method  = LinearBoltzmann::AngleAggregationType::POLAR;
  master_num_grp_subsets = 1;
  master_num_ang_subsets = 1;
  residual_tolerance = 1.0e-6;
  max_iterations = 200;
  gmres_restart_intvl = 30;
  apply_wgdsa = false;
  apply_tgdsa = false;

  wgdsa_solver = nullptr;
  tgdsa_solver = nullptr;

  wgdsa_max_iters = 30;
  tgdsa_max_iters = 30;

  wgdsa_tol = 1.0e-4;
  tgdsa_tol = 1.0e-4;

  wgdsa_verbose = false;
  tgdsa_verbose = false;

  allow_cycles = false;

  log_sweep_events = false;
}

//###################################################################
/**Computes the discrete to moment operator.*/
void LBSGroupset::BuildDiscMomOperator(
  int scatt_order,
  LinearBoltzmann::GeometryType geometry_type)
{
  if (geometry_type == LinearBoltzmann::GeometryType::ONED_SLAB)
  {
    quadrature->BuildDiscreteToMomentOperator(scatt_order,1);
  }
  else if (geometry_type == LinearBoltzmann::GeometryType::TWOD_CARTESIAN)
  {
    quadrature->BuildDiscreteToMomentOperator(scatt_order,2);
  }
  else if (geometry_type == LinearBoltzmann::GeometryType::THREED_CARTESIAN)
  {
    quadrature->BuildDiscreteToMomentOperator(scatt_order,3);
  }
}

//###################################################################
/**Computes the moment to discrete operator.*/
void LBSGroupset::BuildMomDiscOperator(
  int scatt_order,
  LinearBoltzmann::GeometryType geometry_type)
{
  if (geometry_type == LinearBoltzmann::GeometryType::ONED_SLAB)
  {
    quadrature->BuildMomentToDiscreteOperator(scatt_order,1);
  }
  else if (geometry_type == LinearBoltzmann::GeometryType::TWOD_CARTESIAN)
  {
    quadrature->BuildMomentToDiscreteOperator(scatt_order,2);
  }
  else if (geometry_type == LinearBoltzmann::GeometryType::THREED_CARTESIAN)
  {
    quadrature->BuildMomentToDiscreteOperator(scatt_order,3);
  }
}


//###################################################################
/**Constructs the groupset subsets.*/
void LBSGroupset::BuildSubsets()
{
  grp_subsets.clear();
  grp_subset_sizes.clear();

  //=================================== Groupset subsets
  int num_gs_subsets = 1;
  if (master_num_grp_subsets <= groups.size())
    num_gs_subsets = master_num_grp_subsets;

  int gs_subset_size = floor(groups.size()/num_gs_subsets);

  for (int ss=0; ss<num_gs_subsets; ss++)
  {
    int subset_ranki = ss*gs_subset_size;
    int subset_size  = gs_subset_size;

    if (ss == (num_gs_subsets-1))
      subset_size = (int)groups.size() - ss*gs_subset_size;

//    grp_subsets.push_back(GsSubSet(subset_ranki,subset_ranki+subset_size-1));
    grp_subsets.emplace_back(subset_ranki,subset_ranki+subset_size-1);
    grp_subset_sizes.push_back(subset_size);

    chi_log.Log(LOG_0)
    << "Groupset subset " << ss << " "
    << subset_ranki << "->" << subset_ranki+subset_size-1;
  }//for ss

  //=================================== Angle subsets
  ang_subsets_top.clear();
  ang_subsets_bot.clear();
  ang_subset_sizes_top.clear();
  ang_subset_sizes_bot.clear();
  if (quadrature->type == chi_math::AngularQuadratureType::ProductQuadrature)
  {
    auto prodquadrature =
      std::static_pointer_cast<chi_math::ProductQuadrature>(quadrature);
    int num_pol_angls_hemi = (int)prodquadrature->polar_ang.size()/2;
    int num_an_subsets = 1;
    if (master_num_ang_subsets <= num_pol_angls_hemi)
      num_an_subsets = master_num_ang_subsets;

    int an_subset_size = floor(num_pol_angls_hemi/num_an_subsets);

    //==================== Top hemishpere
    for (int ss=0; ss<num_an_subsets; ss++)
    {
      int subset_ranki = ss*an_subset_size;
      int subset_size  = an_subset_size;

      if (ss == (num_an_subsets-1))
        subset_size = num_pol_angls_hemi - ss*an_subset_size;

//      ang_subsets_top.push_back(AngSubSet(subset_ranki,subset_ranki+subset_size-1));
      ang_subsets_top.emplace_back(subset_ranki,subset_ranki+subset_size-1);
      ang_subset_sizes_top.push_back(subset_size);

      if (angleagg_method != LinearBoltzmann::AngleAggregationType::SINGLE)
        chi_log.Log(LOG_0)
          << "Top-hemi Angle subset " << ss << " "
          << subset_ranki << "->" << subset_ranki+subset_size-1;
    }//for ss

    //==================== Bottom hemisphere
    for (int ss=0; ss<num_an_subsets; ss++)
    {
      int subset_ranki = ss*an_subset_size + num_pol_angls_hemi;
      int subset_size  = an_subset_size;

      if (ss == (num_an_subsets-1))
        subset_size = num_pol_angls_hemi - ss*an_subset_size;

//      ang_subsets_bot.push_back(AngSubSet(subset_ranki,subset_ranki+subset_size-1));
      ang_subsets_bot.emplace_back(subset_ranki,subset_ranki+subset_size-1);
      ang_subset_sizes_bot.push_back(subset_size);

      if (angleagg_method != LinearBoltzmann::AngleAggregationType::SINGLE)
        chi_log.Log(LOG_0)
          << "Bot-hemi Angle subset " << ss << " "
          << subset_ranki << "->" << subset_ranki+subset_size-1;
    }//for ss
  }

}

//###################################################################
/**Constructs the groupset subsets.*/
void LBSGroupset::PrintSweepInfoFile(size_t ev_tag, const std::string& file_name)
{
  if (not log_sweep_events) return;

  std::ofstream ofile;
  ofile.open(file_name,std::ofstream::out);

  ofile
    << "Groupset Sweep information "
    << "location " << chi_mpi.location_id << "\n";


  //======================================== Print all anglesets
  for (int q=0; q<angle_agg.angle_set_groups.size(); ++q)
  {
    ofile << "Angle-set group " << q << ":\n";
    auto& ang_set_grp = angle_agg.angle_set_groups[q];
    int num_ang_sets_per_grp = (int)ang_set_grp.angle_sets.size();
    for (int as=0; as<num_ang_sets_per_grp; ++as)
    {
      auto ang_set = ang_set_grp.angle_sets[as];

      int ang_set_num = as + q*num_ang_sets_per_grp;

      ofile << "  Angle-set " << ang_set_num << " angles [# varphi theta]:\n";

      for (auto& ang_num : ang_set->angles)
      {
        const auto& angle = quadrature->abscissae[ang_num];

        ofile
          << "    " << ang_num
          << " " << angle.phi
          << " " << angle.theta << "\n";
      }
    }
  }

  //======================================== Print event history
  ofile << chi_log.PrintEventHistory(ev_tag);

  ofile.close();
}

