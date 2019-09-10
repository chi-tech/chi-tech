#include "lbs_groupset.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>

//###################################################################
/**Computes the discrete to moment operator.*/
void LBSGroupset::BuildDiscMomOperator(int scatt_order)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int num_angles = quadrature->abscissae.size();
  int num_moms = 0;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          double w = quadrature->weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//line mesh
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          double w = quadrature->weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//2d meshers
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          double w = quadrature->weights[n];
          cur_mom.push_back(value*w);
        }

        d2m_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//extruder


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
void LBSGroupset::BuildMomDiscOperator(int scatt_order)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int num_angles = quadrature->abscissae.size();
  int num_moms = 0;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D Slab
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = ((2.0*ell+1.0)/2.0)*
                         chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//line mesh

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = ((2.0*ell+1.0)/2.0/M_PI)*
                         chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//extruder

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
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
          chi_math::QuadraturePointPhiTheta* cur_angle = quadrature->abscissae[n];
          double value = ((2.0*ell+1.0)/4.0/M_PI)*
                         chi_math::Ylm(ell,m,
                                                       cur_angle->phi,
                                                       cur_angle->theta);
          cur_mom.push_back(value);
        }

        m2d_op.push_back(cur_mom);
      }//for m
    }//for ell
  }//extruder

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
  //=================================== Groupset subsets
  int num_gs_subsets = 1;
  if (master_num_grp_subsets < groups.size())
    num_gs_subsets = master_num_grp_subsets;

  int gs_subset_size = ceil(groups.size()/num_gs_subsets);
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
  int num_pol_angls_hemi = quadrature->polar_ang.size()/2;
  int num_an_subsets = 1;
  if (master_num_ang_subsets < num_pol_angls_hemi)
    num_an_subsets = master_num_ang_subsets;

  int an_subset_size = ceil(num_pol_angls_hemi/num_an_subsets);

  //==================== Top hemishpere
  for (int ss=0; ss<num_an_subsets; ss++)
  {
    int subset_ranki = ss*an_subset_size;
    int subset_size  = an_subset_size;

    if (ss == (num_an_subsets-1))
      subset_size = num_pol_angls_hemi - ss*an_subset_size;

    ang_subsets_top.push_back(AngSubSet(subset_ranki,subset_ranki+subset_size-1));
    ang_subset_sizes_top.push_back(subset_size);

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

    chi_log.Log(LOG_0)
      << "Bot-hemi Angle subset " << ss << " "
      << subset_ranki << "->" << subset_ranki+subset_size-1;
  }//for ss
}