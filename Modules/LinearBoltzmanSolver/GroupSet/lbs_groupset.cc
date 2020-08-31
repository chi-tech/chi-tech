#include "lbs_groupset.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>
#include <ChiMesh/VolumeMesher/Predefined3D/volmesher_predefined3d.h>

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
  iterative_method = NPT_GMRES;
  angleagg_method  = LinearBoltzman::AngleAggregationType::POLAR;
  angle_agg = new AngleAgg;
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

  latest_convergence_metric = 1.0;
}

//###################################################################
/**Computes the discrete to moment operator.*/
void LBSGroupset::BuildDiscMomOperator(
  int scatt_order,
  LinearBoltzman::GeometryType geometry_type)
{
  int num_angles = quadrature->abscissae.size();
  int num_moms = 0;

  d2m_op.clear();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (geometry_type == LinearBoltzman::GeometryType::ONED_SLAB)
  {
    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=0; m<=0; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = quadrature->abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          double w = quadrature->weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//line mesh
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D and 3D
  else if ( (geometry_type == LinearBoltzman::GeometryType::TWOD_CARTESIAN) or
            (geometry_type == LinearBoltzman::GeometryType::THREED_CARTESIAN) )
  {
    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=-ell; m<=ell; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = quadrature->abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          double w = quadrature->weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//extruder
  else
  {
    chi_log.Log(LOG_ALLERROR) << "Unsupported geometry type encountered in "
                                 "call to LBSGroupset::BuildDiscMomOperator.";
    exit(EXIT_FAILURE);
  }


  std::stringstream outs;
  outs
    << "\nQuadrature d2m operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << d2m_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi_log.Log(LOG_0VERBOSE_1) << outs.str();
}

//###################################################################
/**Computes the moment to discrete operator.*/
void LBSGroupset::BuildMomDiscOperator(
  int scatt_order,
  LinearBoltzman::GeometryType geometry_type)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int num_angles = quadrature->abscissae.size();
  int num_moms = 0;

  m2d_op.clear();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (geometry_type == LinearBoltzman::GeometryType::ONED_SLAB)
  {
    int mc=-1;
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=0; m<=0; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = quadrature->abscissae[n];
          double value = ((2.0*ell+1.0)/2.0)*
                         chi_math::Ylm(ell,m,
                                       cur_angle.phi,
                                       cur_angle.theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//line mesh
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D and 3D
  else if ( (geometry_type == LinearBoltzman::GeometryType::TWOD_CARTESIAN) or
            (geometry_type == LinearBoltzman::GeometryType::THREED_CARTESIAN) )
  {
    int mc=-1;
    for (int ell=0; ell<=scatt_order; ell++)
    {
      for (int m=-ell; m<=ell; m++)
      {
        std::vector<double> cur_mom; mc++;
        num_moms++;

        for (int n=0; n<num_angles; n++)
        {
          const auto& cur_angle = quadrature->abscissae[n];
          double value = ((2.0*ell+1.0)/4.0/M_PI)*
                         chi_math::Ylm(ell,m,
                         cur_angle.phi,
                         cur_angle.theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//extruder
  else
  {
    chi_log.Log(LOG_ALLERROR) << "Unsupported mesh type encountered in call to"
                                 "LBSGroupset::BuildMomDiscOperator.";
    exit(EXIT_FAILURE);
  }

  std::stringstream outs;

  outs
    << "\nQuadrature m2d operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
        << std::setw(15) << std::left << std::fixed
        << std::setprecision(10) << m2d_op[m][n] << " ";
    }
    outs << "\n";
  }

  chi_log.Log(LOG_0VERBOSE_1) << outs.str();
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
      subset_size = groups.size() - ss*gs_subset_size;

    grp_subsets.push_back(GsSubSet(subset_ranki,subset_ranki+subset_size-1));
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
    int num_pol_angls_hemi = prodquadrature->polar_ang.size()/2;
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

      ang_subsets_top.push_back(AngSubSet(subset_ranki,subset_ranki+subset_size-1));
      ang_subset_sizes_top.push_back(subset_size);

      if (angleagg_method != LinearBoltzman::AngleAggregationType::SINGLE)
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

      ang_subsets_bot.push_back(AngSubSet(subset_ranki,subset_ranki+subset_size-1));
      ang_subset_sizes_bot.push_back(subset_size);

      if (angleagg_method != LinearBoltzman::AngleAggregationType::SINGLE)
        chi_log.Log(LOG_0)
          << "Bot-hemi Angle subset " << ss << " "
          << subset_ranki << "->" << subset_ranki+subset_size-1;
    }//for ss
  }

}

//###################################################################
/**Constructs the groupset subsets.*/
void LBSGroupset::PrintSweepInfoFile(size_t ev_tag, std::string file_name)
{
  if (not log_sweep_events) return;

  std::ofstream ofile;
  ofile.open(file_name,std::ofstream::out);

  ofile
    << "Groupset Sweep information "
    << "location " << chi_mpi.location_id << "\n";


  //======================================== Print all anglesets
  for (int q=0; q<angle_agg->angle_set_groups.size(); ++q)
  {
    ofile << "Angle-set group " << q << ":\n";
    auto ang_set_grp = angle_agg->angle_set_groups[q];
    size_t num_ang_sets_per_grp = ang_set_grp->angle_sets.size();
    for (int as=0; as<num_ang_sets_per_grp; ++as)
    {
      auto ang_set = ang_set_grp->angle_sets[as];

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