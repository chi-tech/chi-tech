#include "lbs_linear_boltzmann_solver.h"



//###################################################################
/**Constructor for LBS*/
LinearBoltzmann::Solver::Solver(const std::string& in_text_name) :
  chi_physics::Solver(in_text_name)
{
  boundary_types.resize(6,
    std::pair<BoundaryType,int>(LinearBoltzmann::BoundaryType::VACUUM, -1));
}