#include "mg_diffusion_solver.h"
#include "tools/tools.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

//#include "ChiMesh/MeshHandler/chi_meshhandler.h"
//#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//========================================================== Execute
void mg_diffusion::Solver::Execute()
{
  chi::log.Log() << "\nExecuting CFEM Multigroup Diffusion solver";

  //============================================= Create Krylov Solver
  // setup KSP once for all
  petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
      A.front(),        //Matrix
      TextName(),      //Solver name
      KSPCG,           //Solver type
      PCGAMG,          //Preconditioner type
      basic_options_("residual_tolerance").FloatValue(),  //Relative residual tolerance
      basic_options_("max_inner_iters").IntegerValue()    //Max # of inner iterations
    );

  KSPSetApplicationContext(petsc_solver.ksp, (void*)&my_app_context);
  KSPMonitorCancel(petsc_solver.ksp);
  KSPMonitorSet(petsc_solver.ksp, &mg_diffusion::MGKSPMonitor,
                nullptr, nullptr);

  int64_t iverbose = basic_options_("verbose_level").IntegerValue();
  my_app_context.verbose = iverbose > 1 ? PETSC_TRUE : PETSC_FALSE;
//  if (my_app_context.verbose == PETSC_TRUE)
//    cout << "--context TRUE" << endl;
//  if (my_app_context.verbose == PETSC_FALSE)
//    cout << "--context FALSE" << endl;
//  std::cout << "STOP" << std::endl; std::cin.get();

  // shortcuts
  // unsigned int lfg = mg_diffusion::Solver::last_fast_group;
  // unsigned int ng = mg_diffusion::Solver::num_groups;

  //============================================= Solve fast groups:
  for (unsigned int g=0; g < last_fast_group; ++g)
  {
    mg_diffusion::Solver::Assemble_RHS(g, iverbose);
    mg_diffusion::Solver::SolveOneGroupProblem(g, iverbose);
  }

  //============================================= Solve thermal groups:
  unsigned int thermal_iteration = 0;
  // max # of thermal iterations
  int64_t max_thermal_iters = basic_options_("max_thermal_iters").IntegerValue();
  // max thermal error between two successive iterates
  double  thermal_tol = basic_options_("thermal_flux_tolerance").FloatValue();
  // computed error
  double thermal_error_all;
  double thermal_error_g;

  do
  {
    thermal_error_all = 0.0;
    for (unsigned int g=last_fast_group; g < num_groups; ++g)
    {
      // conpute rhs src
      mg_diffusion::Solver::Assemble_RHS(g, iverbose);
      // copy solution
      VecCopy(x[g], x_old[g]);
      // solve group g for new solution
      mg_diffusion::Solver::SolveOneGroupProblem(g, iverbose);
      // compute L2 norm of thermal error for current g (requires one more copy)
      VecCopy(x[g], thermal_dphi);
      VecAXPY(thermal_dphi, -1.0, x_old[g]);
      VecNorm(thermal_dphi,NORM_2,&thermal_error_g);
      thermal_error_all = std::max(thermal_error_all,thermal_error_g);
    }
    // perform two-grid
    if (do_two_grid)
    {
      mg_diffusion::Solver::Assemble_RHS_TwoGrid(iverbose);
      mg_diffusion::Solver::SolveOneGroupProblem(num_groups, iverbose);
      mg_diffusion::Solver::Update_Flux_With_TwoGrid(iverbose);
    }

    if (iverbose > 0)
      chi::log.Log() << " --thermal iteration = " << std::setw(5)  << std::right << thermal_iteration
                << ", Error=" << std::setw(11) << std::right << std::scientific << std::setprecision(7)
                << thermal_error_all << std::endl;

    ++thermal_iteration;
  }
  while ( (thermal_error_all > thermal_tol) &&
          (thermal_iteration < max_thermal_iters) );

  if (iverbose > 0)
  {
    if (thermal_error_all < thermal_tol)
      std::cout << "\nThermal iterations converged for fixed-source problem" << std::endl;
    else
      std::cout << "\nThermal iterations NOT converged for fixed-source problem" << std::endl;
  }

  UpdateFieldFunctions();
  chi::log.Log() << "Done solving multi-group diffusion";

}