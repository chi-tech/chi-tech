#include "mg_diffusion_solver.h"
#include "tools/tools.h"

#include "chi_runtime.h"
#include "chi_log.h"

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
      basic_options("residual_tolerance").FloatValue(),  //Relative residual tolerance
      basic_options("max_inner_iters").IntegerValue()    //Max # of inner iterations
    );

  KSPSetApplicationContext(petsc_solver.ksp, (void*)&my_app_context);
  KSPMonitorCancel(petsc_solver.ksp);
  KSPMonitorSet(petsc_solver.ksp, &mg_diffusion::MGKSPMonitor,
                nullptr, nullptr);

  int64_t iverbose = basic_options("verbose_level").IntegerValue();
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
  int64_t max_thermal_iters = basic_options("max_thermal_iters").IntegerValue();
  // max thermal error between two successive iterates
  double  thermal_tol = basic_options("thermal_flux_error").FloatValue();
  // computed error
  double thermal_error = 0.0;

  do
  {
    thermal_error = 0.0;
    for (unsigned int g=last_fast_group; g < num_groups; ++g)
    {
      mg_diffusion::Solver::Assemble_RHS(g, iverbose);
      mg_diffusion::Solver::SolveOneGroupProblem(g, iverbose);

      // copy solution
      VecCopy(x[g], thermal_dphi);
      // compute the difference
      VecAXPY(thermal_dphi, -1.0, x_old[g]);
      // compute the L2 norm
      VecNorm(thermal_dphi,NORM_2,&thermal_error);
      // copy solution
      VecCopy(x[g], x_old[g]);
    }
    // two-grid
    if (do_two_grid)
    {
      mg_diffusion::Solver::Assemble_RHS_TwoGrid(iverbose);
      mg_diffusion::Solver::SolveOneGroupProblem(num_groups, iverbose);
    }
    if ( (iverbose > 0) && (num_groups != last_fast_group) )
      std::cout << " --thermal iteration = " << std::setw(6)  << std::right << thermal_iteration
                << ", Error=" << std::setw(11) << std::right << std::scientific << std::setprecision(5)
                << thermal_error << std::endl;

    ++thermal_iteration;
  }
  while ( (num_groups != last_fast_group)
          && (thermal_error > thermal_tol)
          && (thermal_iteration < max_thermal_iters) );

  if ( (iverbose > 0) && (num_groups != last_fast_group) )
  {
    if (thermal_error < thermal_tol)
      std::cout << "Thermal iterations converged for fixed-source problem" << std::endl;
    else
      std::cout << "Thermal iterations NOT converged for fixed-source problem" << std::endl;
  }

  chi::log.Log() << "Done solving multi-group diffusion";

}