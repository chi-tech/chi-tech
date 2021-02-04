#include "lbs_linear_boltzmann_solver.h"



//###################################################################
/**Constructor for LBS*/
LinearBoltzmann::Solver::Solver()
{
  //============================================= Default options
  max_cell_dof_count = 0;

  discretization = nullptr;

  boundary_types.resize(6,
    std::pair<BoundaryType,int>(LinearBoltzmann::BoundaryType::VACUUM, -1));
}