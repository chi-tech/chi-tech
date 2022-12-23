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

  //============================================= Solve fast groups:
  for (unsigned int g=0; g<mg_diffusion::Solver::last_fast_group; ++g)
  {
    mg_diffusion::Solver::Assemble_RHS(g);
    mg_diffusion::Solver::SolveOneGroupProblem(g);
  }
  chi::log.Log() << "Done solving";

}