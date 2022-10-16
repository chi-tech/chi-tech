#include "fv.h"

#include "chi_log.h"

#include "chi_runtime.h"
#include "ChiTimer/chi_timer.h"

//###################################################################
/**Only constructor for this method.*/
chi_math::SpatialDiscretization_FV::
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr& in_grid,
                           chi_math::CoordinateSystemType in_cs_type)
  : chi_math::SpatialDiscretization(in_grid, in_cs_type, SDMType::FINITE_VOLUME)
{
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Creating Finite Volume spatial discretizaiton.";

  PreComputeCellSDValues();
  PreComputeNeighborCellSDValues();
  OrderNodes();
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Done creating Finite Volume spatial discretizaiton.";
}