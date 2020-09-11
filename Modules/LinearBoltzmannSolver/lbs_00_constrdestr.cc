#include "lbs_linear_boltzmann_solver.h"



//###################################################################
/**Constructor for NPT*/
LinearBoltzmann::Solver::Solver()
{
  //============================================= Default options
  options.scattering_order  = 1;
  options.partition_method = PARTITION_METHOD_SERIAL;

  max_cell_dof_count = 0;

  discretization = nullptr;

  boundary_types.resize(6,
    std::pair<BoundaryType,int>(LinearBoltzmann::BoundaryType::VACUUM, -1));
}