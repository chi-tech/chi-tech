#include "gs_linear_solver.h"

#include "gs_convergence_test.h"

#include "LinearBoltzmannSolvers/LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"
#include "ChiMath/LinearSolver/linear_matrix_action_Ax.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <petscksp.h>
#include <memory>
#include <iomanip>

#define sc_double static_cast<double>
#define sc_int64_t static_cast<int64_t>

#define GetGSContextPtr(x) \
        std::dynamic_pointer_cast<GSContext<Mat,Vec,KSP>>(x)

namespace lbs
{

template<>
void GSLinearSolver<Mat, Vec, KSP>::PreSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PreSetupCallback();
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetSolverContext()
{
  KSPSetApplicationContext(solver_, &(*context_ptr_));
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetConvergenceTest()
{
  KSPSetConvergenceTest(solver_, &GSConvergenceTest, nullptr, nullptr);
}


template<>
void GSLinearSolver<Mat, Vec, KSP>::SetSystemSize()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);
  const auto sizes = gs_context_ptr->SystemSize();

  num_local_dofs_ = sizes.first;
  num_globl_dofs_ = sizes.second;
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetSystem()
{
  x_ = chi_math::PETScUtils::CreateVector(sc_int64_t(num_local_dofs_),
                                          sc_int64_t(num_globl_dofs_));

  VecSet(x_,0.0);
  VecDuplicate(x_,&b_);

  //============================================= Create the matrix-shell
  MatCreateShell(PETSC_COMM_WORLD,sc_int64_t(num_local_dofs_),
                                  sc_int64_t(num_local_dofs_),
                                  sc_int64_t(num_globl_dofs_),
                                  sc_int64_t(num_globl_dofs_),
                                  &(*context_ptr_),&A_);

  //============================================= Set the action-operator
  MatShellSetOperation(A_, MATOP_MULT,
                     (void (*)()) chi_math::LinearSolverMatrixAction<Mat, Vec>);

  //============================================= Set solver operators
  KSPSetOperators(solver_, A_, A_);
  KSPSetUp(solver_);
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetPreconditioner()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->SetPreconditioner(solver_);
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::PostSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PostSetupCallback();
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::PreSolveCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PreSolveCallback();
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetRHS()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset   = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  if (gs_context_ptr->log_info_)
    chi::log.Log() << chi::program_timer.GetTimeString() << " Computing b";

  //SetSource for RHS
  saved_q_moments_local_ = lbs_solver.QMomentsLocal();
  const int scope = gs_context_ptr->rhs_src_scope_;
  gs_context_ptr->set_source_function_(groupset,
                                       lbs_solver.QMomentsLocal(),
                                       lbs_solver.PhiOldLocal(),
                                       scope);

  //=================================================== Apply transport operator
  gs_context_ptr->ApplyInverseTransportOperator(scope);

  //=================================================== Assemble PETSc vector
  lbs_solver.
    SetGSPETScVecFromPrimarySTLvector(groupset,
                                      b_,
                                      lbs_solver.PhiNewLocal(),
                                      gs_context_ptr->with_delayed_psi_);

  //============================================= Compute RHS norm
  VecNorm(b_, NORM_2, &context_ptr_->rhs_norm_);

  //============================================= Compute precondition RHS norm
  PC pc;
  KSPGetPC(solver_, &pc);
  Vec temp_vec;
  VecDuplicate(b_, &temp_vec);
  PCApply(pc, b_, temp_vec);
  VecNorm(temp_vec, NORM_2, &context_ptr_->rhs_preconditioned_norm_);
  VecDestroy(&temp_vec);
}

/**Sets the initial guess for a gs solver. If the initial guess's norm
 * is large enough the initial guess will be used, otherwise it is assumed
 * zero.*/
template<> void GSLinearSolver<Mat, Vec, KSP>::SetInitialGuess()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset   = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  lbs_solver.
  SetGSPETScVecFromPrimarySTLvector(groupset, x_, lbs_solver.PhiOldLocal());

  double init_guess_norm = 0.0;
  VecNorm(x_,NORM_2,&init_guess_norm);

  if (init_guess_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(solver_, PETSC_TRUE);
    if (gs_context_ptr->log_info_)
      chi::log.Log() << "Using phi_old as initial guess.";
  }
}

/**For this callback we simply restore the q_moments_local vector.*/
template<> void GSLinearSolver<Mat, Vec, KSP>::PostSolveCallback()
{
  //============================================= Get convergence reason
  KSPConvergedReason reason;
  KSPGetConvergedReason(solver_,&reason);
  if (reason != KSP_CONVERGED_RTOL)
    chi::log.Log0Warning()
      << "Krylov solver failed. "
      << "Reason: " << chi_physics::GetPETScConvergedReasonstring(reason);

  //============================================= Copy x to local solution
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset   = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  auto& phi_new_local = lbs_solver.PhiNewLocal();
  auto& phi_old_local = lbs_solver.PhiOldLocal();

  const bool with_delayed_psi = gs_context_ptr->with_delayed_psi_;

  lbs_solver.
    SetPrimarySTLvectorFromGSPETScVec(groupset, x_,
                                      phi_new_local, with_delayed_psi);
  lbs_solver.
    SetPrimarySTLvectorFromGSPETScVec(groupset, x_,
                                      phi_old_local, with_delayed_psi);

  //============================================= Restore saved q_moms
  lbs_solver.QMomentsLocal() = saved_q_moments_local_;

  //============================================= Context specific callback
  gs_context_ptr->PostSolveCallback();
}
}//namespace lbs

