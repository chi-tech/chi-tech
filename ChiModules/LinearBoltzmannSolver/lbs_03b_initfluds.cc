#include "lbs_linear_boltzmann_solver.h"
#include "ChiConsole/chi_console.h"

#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include <iomanip>

extern ChiConsole&  chi_console;
extern ChiLog&     chi_log;
extern ChiMPI&      chi_mpi;

#include "ChiTimer/chi_timer.h"

extern ChiTimer chi_program_timer;

//###################################################################
/**Initializes fluds data structures.*/
void LinearBoltzmann::Solver::InitFluxDataStructures(LBSGroupset& groupset)
{
  //================================================== Build angular flux unknown
  groupset.psi_uk_man.unknowns.clear();
  size_t num_angles = groupset.quadrature->abscissae.size();
  size_t num_groups = groupset.groups.size();
  auto& grpset_psi_uk_man = groupset.psi_uk_man;

  for (unsigned int n=0; n<num_angles; ++n)
    grpset_psi_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);

  if (options.save_angular_flux)
  {
    size_t num_ang_unknowns = discretization->GetNumLocalDOFs(grid,
                                                              grpset_psi_uk_man);
    groupset.psi_to_be_saved = true;
    groupset.num_psi_unknowns_local = num_ang_unknowns;
    groupset.psi_new_local.assign(num_ang_unknowns,0.0);
  }

  //================================================== Angle Aggregation
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher& mesher = *handler->volume_mesher;

  if ( options.geometry_type == GeometryType::ONED_SLAB or
       options.geometry_type == GeometryType::TWOD_CARTESIAN or
       (typeid(mesher) == typeid(chi_mesh::VolumeMesherExtruder)))
  {
    if      (groupset.angleagg_method == AngleAggregationType::SINGLE)
    {
      InitAngleAggSingle(groupset);
    }
    else if (groupset.angleagg_method == AngleAggregationType::POLAR)
    {
      InitAngleAggPolar(groupset);
    }
  }
  else
  {
    InitAngleAggSingle(groupset);
  }

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Initialized Angle Aggregation.   "
    << "         Process memory = "
    << std::setprecision(3) << chi_console.GetMemoryUsageInMB()
    << " MB.";


  MPI_Barrier(MPI_COMM_WORLD);
}