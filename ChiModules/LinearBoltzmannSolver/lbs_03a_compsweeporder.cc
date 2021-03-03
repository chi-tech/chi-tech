#include "lbs_linear_boltzmann_solver.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "ChiMath/Quadratures/product_quadrature.h"

#include "chi_mpi.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;
extern ChiTimer chi_program_timer;

typedef chi_mesh::sweep_management::AngleSet TAngleSet;
typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;

#include <iomanip>
#include "ChiConsole/chi_console.h"

extern ChiConsole&  chi_console;

//###################################################################
/**Initializes the sweep ordering for the given groupset.*/
void LinearBoltzmann::Solver::ComputeSweepOrderings(LBSGroupset& groupset) const
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Computing Sweep ordering.\n";

  //============================================= Clear sweep ordering
  groupset.sweep_orderings.clear();
  groupset.sweep_orderings.shrink_to_fit();

  auto mesh_handler = chi_mesh::GetCurrentHandler();
  auto mesher = mesh_handler->volume_mesher;

  //============================================= Check possibility of cycles
  if (mesher->options.partition_type ==
      chi_mesh::VolumeMesher::PartitionType::PARMETIS and
      not groupset.allow_cycles)
  {
    chi_log.Log(LOG_ALLERROR)
      << "When using PARMETIS type partitioning then groupset iterative method"
         " must be NPT_CLASSICRICHARDSON_CYCLES or NPT_GMRES_CYCLES";
    exit(EXIT_FAILURE);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single angle aggr.
  if (groupset.angleagg_method == LinearBoltzmann::AngleAggregationType::SINGLE)
  {
    for (auto& angle : groupset.quadrature->abscissae)
    {
      auto new_swp_order =
        chi_mesh::sweep_management::
        CreateSweepOrder(angle.theta,
                         angle.phi,
                         this->grid,
                         groupset.allow_cycles);
      groupset.sweep_orderings.push_back(new_swp_order);
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D MESHES
  else if (options.geometry_type == GeometryType::ONED_SLAB)
  {
    if (groupset.quadrature->type == chi_math::AngularQuadratureType::ProductQuadrature)
    {
      auto product_quadrature =
        std::static_pointer_cast<chi_math::ProductQuadrature>(groupset.quadrature);

      int num_azi = product_quadrature->azimu_ang.size();
      int num_pol = product_quadrature->polar_ang.size();
      int pa      = num_pol/2;

      if (num_azi != 1)
      {
        chi_log.Log(LOG_0)
          << "Incompatible number of azimuthal angles in quadrature set "
          << "for a 1D simulation.";
        exit(EXIT_FAILURE);
      }

      auto new_swp_order =
        chi_mesh::sweep_management::
        CreateSweepOrder(product_quadrature->polar_ang[0],
                         product_quadrature->azimu_ang[0],
                         this->grid,
                         groupset.allow_cycles);
      groupset.sweep_orderings.push_back(new_swp_order);

      new_swp_order =
        chi_mesh::sweep_management::
        CreateSweepOrder(product_quadrature->polar_ang[pa],
                         product_quadrature->azimu_ang[0],
                         this->grid,
                         groupset.allow_cycles);
      groupset.sweep_orderings.push_back(new_swp_order);
    }

  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D 3D MESHES
  else if ( options.geometry_type == GeometryType::TWOD_CARTESIAN or
            (options.geometry_type == GeometryType::THREED_CARTESIAN and
              typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))/*
            (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder)) or
            (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D)) or
            (typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined3D))*/)
  {
    if (groupset.quadrature->type == chi_math::AngularQuadratureType::ProductQuadrature)
    {
      auto product_quadrature =
        std::static_pointer_cast<chi_math::ProductQuadrature>(groupset.quadrature);


      int num_azi = product_quadrature->azimu_ang.size();
      int num_pol = product_quadrature->polar_ang.size();

      if (num_azi < 4)
      {
        chi_log.Log(LOG_0)
          << "Incompatible number of azimuthal angles in quadrature set "
          << "for a 2D or 3D simulation.";
        exit(EXIT_FAILURE);
      }
      if (num_pol < 2)
      {
        chi_log.Log(LOG_0)
          << "Incompatible number of polar angles in quadrature set "
          << "for a 2D or 3D simulation.";
        exit(EXIT_FAILURE);
      }

      //============================================= Create sweep ordering
      //                                              per azimuthal angle
      //                                              per hemisphere
      int pa = num_pol/2;

      //=========================================== TOP HEMISPHERE
      for (int i=0; i<num_azi; i++)
      {
        auto new_swp_order =
          chi_mesh::sweep_management::
          CreateSweepOrder(product_quadrature->polar_ang[pa-1],
                           product_quadrature->azimu_ang[i],
                           this->grid,
                           groupset.allow_cycles);
        groupset.sweep_orderings.push_back(new_swp_order);
      }
      //=========================================== BOTTOM HEMISPHERE
      for (int i=0; i<num_azi; i++)
      {
        auto new_swp_order =
          chi_mesh::sweep_management::
          CreateSweepOrder(product_quadrature->polar_ang[pa],
                           product_quadrature->azimu_ang[i],
                           this->grid,
                           groupset.allow_cycles);
        groupset.sweep_orderings.push_back(new_swp_order);
      }
    }//if product quadrature
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "The simulation is not using \"LBSGroupset.ANGLE_AGG_SINGLE\", "
           "and therefore only certain angular quadrature types are supported. "
           "i.e., for now just AngularQuadratureType::ProductQuadrature.";
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "The simulation is not using \"LBSGroupset.ANGLE_AGG_SINGLE\", "
         "and therefore only certain geometry types are supported. i.e., "
         "GeometryType::ONED_SLAB, GeometryType::TWOD_CARTESIAN, "
         "GeometryType::THREED_CARTESIAN.";
    exit(EXIT_FAILURE);
  }


  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Done computing sweep orderings.           Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";

}
