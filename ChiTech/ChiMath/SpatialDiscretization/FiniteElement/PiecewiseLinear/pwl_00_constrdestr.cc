#include "pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWLD::
  SpatialDiscretization_PWLD(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder,
                             chi_math::CoordinateSystemType in_cs_type) :
  SpatialDiscretization_FE(0, in_grid, in_cs_type,
                           SDMType::PIECEWISE_LINEAR_DISCONTINUOUS,
                           setup_flags),
  line_quad_order_arbitrary(qorder),
  tri_quad_order_arbitrary(qorder),
  quad_quad_order_arbitrary(qorder),
  tet_quad_order_arbitrary(qorder),
  hex_quad_order_arbitrary(qorder)
{
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Creating Piecewise Linear Discontinuous "
                   "Finite Element spatial discretizaiton.";

  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Communicating partition neighbors.";
  ref_grid->CommunicatePartitionNeighborCells(neighbor_cells);

  if (setup_flags != chi_math::finite_element::NO_FLAGS_SET)
  {
    PreComputeCellSDValues();
    PreComputeNeighborCellSDValues();
  }

  OrderNodes();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done creating Piecewise Linear Discontinuous "
                   "Finite Element spatial discretizaiton.";
}

