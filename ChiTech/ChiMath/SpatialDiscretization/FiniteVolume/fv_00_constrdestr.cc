#include "fv.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Only constructor for this method.*/
SpatialDiscretization_FV::
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr& in_grid)
  : SpatialDiscretization(0, in_grid, SDMType::FINITE_VOLUME)
{
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Creating Finite Volume spatial discretizaiton.";
  mapping_initialized = false;
  PreComputeCellSDValues();
  PreComputeNeighborCellSDValues();
  OrderNodes();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done creating Finite Volume spatial discretizaiton.";
}