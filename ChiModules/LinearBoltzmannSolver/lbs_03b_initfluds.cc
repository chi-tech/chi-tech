#include "lbs_linear_boltzmann_solver.h"
#include "ChiConsole/chi_console.h"

#include <ChiMesh/SweepUtilities/FLUDS/FLUDS.h>
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>

#include <chi_log.h>
#include <chi_mpi.h>

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
  //================================================== Angle Aggregation
  chi_mesh::MeshHandler* handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher& mesher = *handler->volume_mesher;

  if ((typeid(mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D)) or
      (typeid(mesher) == typeid(chi_mesh::VolumeMesherPredefined2D)) or
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