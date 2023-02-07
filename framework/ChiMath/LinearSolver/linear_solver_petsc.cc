#include "linear_solver.h"

#include <petscksp.h>

namespace chi_math
{
template<> void LinearSolver<Mat, Vec, KSP>::PreSetupCallback();

template<> void LinearSolver<Mat, Vec, KSP>::SetOptions();

template<> void LinearSolver<Mat, Vec, KSP>::SetSolverContext();

template<> void LinearSolver<Mat, Vec, KSP>::SetConvergenceTest();

template<> void LinearSolver<Mat, Vec, KSP>::SetMonitor();

template<> void LinearSolver<Mat, Vec, KSP>::SetPreconditioner();

template<> void LinearSolver<Mat, Vec, KSP>::PostSetupCallback();

template<> void LinearSolver<Mat, Vec, KSP>::PreSolveCallback();

template<> void LinearSolver<Mat, Vec, KSP>::PostSolveCallback();

template<>
void LinearSolver<Mat, Vec, KSP>::Setup()
{
  this->PreSetupCallback();

  KSPCreate(PETSC_COMM_WORLD, &solver_);
  KSPSetType(solver_, iterative_method_.c_str());
  KSPSetTolerances(solver_,
                   tolerance_options_.residual_relative_,
                   tolerance_options_.residual_absolute_,
                   tolerance_options_.residual_divergence_,
                   tolerance_options_.maximum_iterations_);

  if (iterative_method_ == "gmres")
  {
    KSPGMRESSetRestart(solver_, tolerance_options_.gmres_restart_interval);
    KSPGMRESSetBreakdownTolerance(solver_,
                                  tolerance_options_.gmres_breakdown_tolerance);
  }

  KSPSetInitialGuessNonzero(solver_, PETSC_FALSE);

  this->SetOptions();

  this->SetSolverContext();
  this->SetConvergenceTest();
  this->SetMonitor();

  this->SetSystemSize();
  this->SetSystem();

  this->SetPreconditioner();

  this->PostSetupCallback();
}

template<>
void LinearSolver<Mat, Vec, KSP>::PreSetupCallback()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::SetOptions()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::SetSolverContext()
{
  context_ptr_ = std::make_shared<LinearSolverContext<Mat,Vec>>();
  KSPSetApplicationContext(solver_, &(*context_ptr_));
}

template<>
void LinearSolver<Mat, Vec, KSP>::SetConvergenceTest()
{
  KSPSetConvergenceTest(solver_, &KSPConvergedDefault, nullptr, nullptr);
}

template<>
void LinearSolver<Mat, Vec, KSP>::SetMonitor()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::SetPreconditioner()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::PostSetupCallback()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::Solve()
{
  this->PreSolveCallback();
  this->SetRHS();
  this->SetInitialGuess();

  KSPSolve(solver_, b_, x_);
  this->PostSolveCallback();
}

template<>
void LinearSolver<Mat, Vec, KSP>::PreSolveCallback()
{}

template<>
void LinearSolver<Mat, Vec, KSP>::PostSolveCallback()
{}


}//namespace chi_math