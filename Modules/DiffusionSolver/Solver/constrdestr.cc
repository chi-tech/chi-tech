#include "diffusion_solver.h"

chi_diffusion::Solver::Solver()
{
  solver_name = std::string("Diffusion Solver");
  verbose = true;
  common_items_initialized = false;
  max_iters = 500;
  residual_tolerance = 1.0e-8;

  gi = 0;
  G = 1;

  fem_method = 0;

  property_map_D = 0;
  property_map_q = 1;
  property_map_sigma = 2;
  material_mode = DIFFUSION_MATERIALS_REGULAR;
  D_field = nullptr;
  q_field = nullptr;
  sigma_field = nullptr;

}

chi_diffusion::Solver::Solver(std::string in_solver_name):Solver()
{
  solver_name = in_solver_name;
}