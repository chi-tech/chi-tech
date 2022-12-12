#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

#include "../lbs_make_subset.h"

#include "ChiMath/Quadratures/angular_product_quadrature.h"

#include <chi_log.h>
#include <chi_mpi.h>

#include <fstream>

//##############################################
/**Groupset constructor.*/
lbs::LBSGroupset::LBSGroupset(int in_id) :
  id(in_id)
{
}

//###################################################################
/**Computes the discrete to moment operator.*/
void lbs::LBSGroupset::BuildDiscMomOperator(
  unsigned int scattering_order,
  lbs::GeometryType geometry_type)
{
  if (geometry_type == lbs::GeometryType::ONED_SLAB ||
      geometry_type == lbs::GeometryType::ONED_CYLINDRICAL ||
      geometry_type == lbs::GeometryType::ONED_SPHERICAL)
  {
    quadrature->BuildDiscreteToMomentOperator(scattering_order,1);
  }
  else if (geometry_type == lbs::GeometryType::TWOD_CARTESIAN ||
           geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
  {
    quadrature->BuildDiscreteToMomentOperator(scattering_order,2);
  }
  else if (geometry_type == lbs::GeometryType::THREED_CARTESIAN)
  {
    quadrature->BuildDiscreteToMomentOperator(scattering_order,3);
  }
}

//###################################################################
/**Computes the moment to discrete operator.*/
void lbs::LBSGroupset::BuildMomDiscOperator(
  unsigned int scattering_order,
  lbs::GeometryType geometry_type)
{
  if (geometry_type == lbs::GeometryType::ONED_SLAB ||
      geometry_type == lbs::GeometryType::ONED_CYLINDRICAL ||
      geometry_type == lbs::GeometryType::ONED_SPHERICAL)
  {
    quadrature->BuildMomentToDiscreteOperator(scattering_order,1);
  }
  else if (geometry_type == lbs::GeometryType::TWOD_CARTESIAN ||
           geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
  {
    quadrature->BuildMomentToDiscreteOperator(scattering_order,2);
  }
  else if (geometry_type == lbs::GeometryType::THREED_CARTESIAN)
  {
    quadrature->BuildMomentToDiscreteOperator(scattering_order,3);
  }
}


//###################################################################
/**Constructs the groupset subsets.*/
void lbs::LBSGroupset::BuildSubsets()
{
  grp_subset_infos = lbs::MakeSubSets(groups.size(), master_num_grp_subsets);
  {
    size_t ss=0;
    for (const auto& info : grp_subset_infos)
    {
      chi::log.Log()
        << "Groupset subset " << ss << " "
      << info.ss_begin << "->" << info.ss_end;
      ++ss;
    }
  }
}

//###################################################################
/**Constructs the groupset subsets.*/
void lbs::LBSGroupset::PrintSweepInfoFile(size_t ev_tag, const std::string& file_name)
{
  if (not log_sweep_events) return;

  std::ofstream ofile;
  ofile.open(file_name,std::ofstream::out);

  ofile
    << "Groupset Sweep information "
    << "location " << chi::mpi.location_id << "\n";


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
  ofile << chi::log.PrintEventHistory(ev_tag);

  ofile.close();
}

