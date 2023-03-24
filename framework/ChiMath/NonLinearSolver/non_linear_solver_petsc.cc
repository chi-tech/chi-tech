#include "non_linear_solver.h"

#include <petscsnes.h>

namespace chi_math
{

template<>
NonLinearSolver<Mat, Vec, SNES>::~NonLinearSolver()
{
  SNESDestroy(&nl_solver_);
  VecDestroy(&x_);
  VecDestroy(&r_);
  MatDestroy(&J_);
}

template<>
void NonLinearSolver<Mat, Vec, SNES>::ApplyToleranceOptions()
{
  SNESSetTolerances(nl_solver_,
                    tolerance_options_.nl_absolute_tol,
                    tolerance_options_.nl_relative_tol,
                    tolerance_options_.nl_solution_tol,
                    tolerance_options_.nl_max_iterations,
                    tolerance_options_.nl_max_r_evaluations);
  SNESSetMaxLinearSolveFailures(nl_solver_,
                                tolerance_options_.l_max_failed_iterations);
  KSP ksp;
  SNESGetKSP(nl_solver_, &ksp);
  KSPSetTolerances(ksp,
                   tolerance_options_.l_relative_tol,
                   tolerance_options_.l_absolute_tol,
                   tolerance_options_.l_divergence_tol,
                   tolerance_options_.l_max_iterations);
  if (linear_method_ == "gmres")
  {
    KSPGMRESSetRestart(ksp, tolerance_options_.l_gmres_restart_interval);
    KSPGMRESSetBreakdownTolerance(ksp,
                                  tolerance_options_.l_gmres_breakdown_tol);
  }
  KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
}
template<>
void NonLinearSolver<Mat, Vec, SNES>::PreSetupCallback()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::SetOptions()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::SetSolverContext()
{
  SNESSetApplicationContext(nl_solver_, &(*context_ptr_));
}

template<>
void NonLinearSolver<Mat, Vec, SNES>::SetConvergenceTest()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::SetMonitor()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::SetPreconditioner()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::PostSetupCallback()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::Setup()
{
  if (IsSystemSet()) return;
  this->PreSetupCallback();

  SNESCreate(PETSC_COMM_WORLD, &nl_solver_);
  SNESSetType(nl_solver_, nl_method_.c_str());

  KSP ksp;
  SNESGetKSP(nl_solver_, &ksp);
  KSPSetType(ksp, linear_method_.c_str());

  this->ApplyToleranceOptions();

  this->SetOptions();

  this->SetSolverContext();
  this->SetOptions();

  this->SetSolverContext();
  this->SetConvergenceTest();
  this->SetMonitor();


  this->SetSystemSize();
  this->SetSystem();

  this->SetFunction();
  this->SetJacobian();

  this->SetPreconditioner();

  this->PostSetupCallback();
  system_set_ = true;
}

template<>
void NonLinearSolver<Mat, Vec, SNES>::PreSolveCallback()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::PostSolveCallback()
{}

template<>
void NonLinearSolver<Mat, Vec, SNES>::Solve()
{
  this->PreSolveCallback();
  this->SetInitialGuess();

  SNESSolve(nl_solver_, nullptr, x_);
  this->PostSolveCallback();
}

}//namespace chi_math
