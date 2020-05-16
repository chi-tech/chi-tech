#include "diffusion_solver.h"

#include <chi_log.h>
extern ChiLog& chi_log;

chi_diffusion::Solver::Solver()
{
  solver_name = std::string("Diffusion Solver");
  verbose_info = true;
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

chi_diffusion::Solver::~Solver()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Cleaning up diffusion solver: " << solver_name;
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
  KSPDestroy(&ksp);

  for (auto ip_cell_view : ip_cell_views)
    delete ip_cell_view;

  for (auto& loc_border_info : ip_locI_bordercell_info)
    for (auto val : loc_border_info)
      delete val;

  for (auto& loc_border_cell : ip_locI_bordercells)
    for (auto val : loc_border_cell)
      delete val;

  for (auto& loc_fe_view : ip_locI_borderfeviews)
    for (auto val : loc_fe_view)
      delete val;

  for (auto& loc_border_ipview : ip_locI_borderipviews)
    for (auto val : loc_border_ipview)
      delete val;

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Done cleaning up diffusion solver: " << solver_name;
}