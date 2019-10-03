#include "lbs_linear_boltzman_solver.h"



//###################################################################
/**Constructor for NPT*/
LinearBoltzman::Solver::Solver()
{
  //============================================= Default options
  options.scattering_order  = 1;
  options.partition_method = PARTITION_METHOD_SERIAL;

  max_cell_dof_count = 0;

  discretization = nullptr;

  boundary_types.resize(6,std::pair<int,int>(VACUUM,-1));
}