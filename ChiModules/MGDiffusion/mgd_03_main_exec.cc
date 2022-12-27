#include "mg_diffusion_solver.h"

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
      basic_options("max_iters").IntegerValue()          //Max iterations
    );

  // shortcuts
  unsigned int lfg = mg_diffusion::Solver::last_fast_group;
  unsigned int ng = mg_diffusion::Solver::num_groups;
  int i_=0;
  cout << i_<<std::endl; i_++;

  //============================================= Solve fast groups:
  for (unsigned int g=0; g<lfg; ++g)
  {
    mg_diffusion::Solver::Assemble_RHS(g);
    mg_diffusion::Solver::SolveOneGroupProblem(g);
  }

  //============================================= Solve thermal groups:
  unsigned int thermal_iteration = 0;
  double thermal_error = 0.0;

  bool verbose = true;
  int counter=0;
  do
  {
    thermal_error = 0.0;
    for (unsigned int g=lfg; g<ng; ++g)
    {
      mg_diffusion::Solver::Assemble_RHS(g);
      mg_diffusion::Solver::SolveOneGroupProblem(g);

      // copy solution
      VecCopy(x[g], thermal_dphi);
      // compute the difference
      VecAXPY(thermal_dphi, -1.0, x_old[g]);
      // compute the L2 norm
      VecNorm(thermal_dphi,NORM_2,&thermal_error);
      // copy solution
      VecCopy(x[g], x_old[g]);
    }
    if ( verbose && (ng != lfg) )
      std::cout << " --thermal iteration = " << std::setw(6)  << std::right << thermal_iteration
                << "Error=" << std::setw(11) << std::right << std::scientific << std::setprecision(5)
                << thermal_error << std::endl;

    ++thermal_iteration;
  }
  while ((ng!=lfg) && (thermal_error > 1e-6) && (thermal_iteration < 100) );

  if (  verbose && (ng != lfg) )
      std::cout << "Thermal iterations converged for fixed-source problem" << std::endl;


  chi::log.Log() << "Done solving";

}