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
#define WITH_DELAYED_PSI true

namespace lbs
{

template<>
void GSLinearSolver<Mat, Vec, KSP>::PreSetupCallback()
{
  if (log_info)
  {
    std::string method_name;
    switch (groupset_.iterative_method)
    {
      case IterativeMethod::KRYLOV_RICHARDSON:
        method_name = "KRYLOV_RICHARDSON"; break;
      case IterativeMethod::KRYLOV_GMRES:
        method_name = "KRYLOV_GMRES"; break;
      case IterativeMethod::KRYLOV_BICGSTAB:
        method_name = "KRYLOV_BICGSTAB"; break;
      default: method_name = "KRYLOV_GMRES";
    }
    chi::log.Log()
      << "\n\n";
    chi::log.Log()
      << "********** Solving groupset " << groupset_.id
      << " with " << method_name << ".\n\n";
    chi::log.Log()
      << "Quadrature number of angles: "
      << groupset_.quadrature->abscissae.size() << "\n"
      << "Groups " << groupset_.groups.front().id << " "
      << groupset_.groups.back().id << "\n\n";
  }
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetSolverContext()
{
  context_ptr_ = std::make_shared<GSContext<Mat,Vec>>(
    groupset_,
    lbs_solver_,
    set_source_function_,
    lhs_src_scope_,
    rhs_src_scope_);

  KSPSetApplicationContext(solver_, &(*context_ptr_));
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetConvergenceTest()
{
  KSPSetConvergenceTest(solver_, &GSConvergenceTest, nullptr, nullptr);
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetSystemMatrix()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();
  const size_t num_moments      = lbs_solver_.NumMoments();

  const size_t groupset_numgrps = groupset_.groups.size();
  const auto num_delayed_psi_info = groupset_.angle_agg.GetNumDelayedAngularDOFs();
  const size_t local_size = local_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.first;
  const size_t globl_size = globl_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.second;
  const size_t num_angles = groupset_.quadrature->abscissae.size();
  const size_t num_psi_global = globl_node_count *
                                num_angles *
                                groupset_.groups.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info)
  {
    chi::log.Log()
      << "Total number of angular unknowns: "
      << num_psi_global
      << "\n"
      << "Number of lagged angular unknowns: "
      << num_delayed_psi_globl << "("
      << std::setprecision(2)
      << sc_double(num_delayed_psi_globl)*100 / sc_double(num_psi_global)
      << "%)";
  }

  x_ = chi_math::PETScUtils::CreateVector(sc_int64_t(local_size),
                                          sc_int64_t(globl_size));

  VecSet(x_,0.0);
  VecDuplicate(x_,&b_);

  //============================================= Create the matrix-shell
  MatCreateShell(PETSC_COMM_WORLD,sc_int64_t(local_size),
                                  sc_int64_t(local_size),
                                  sc_int64_t(globl_size),
                                  sc_int64_t(globl_size),
                                  &(*context_ptr_),&A_);

  //============================================= Set the action-operator
  MatShellSetOperation(A_, MATOP_MULT,
                     (void (*)()) chi_math::LinearSolverMatrixAction<Mat, Vec>);

  //============================================= Set solver operators
  KSPSetOperators(solver_, A_, A_);
}

template<>
void GSLinearSolver<Mat, Vec, KSP>::SetRHS()
{
  if (log_info)
  {
    chi::log.Log() << chi::program_timer.GetTimeString() << " Computing b";
  }

  //SetSource for RHS
  saved_q_moments_local_ = lbs_solver_.QMomentsLocal();
  set_source_function_(groupset_, lbs_solver_.QMomentsLocal(), rhs_src_scope_);

  //=================================================== Apply transport operator
  auto gs_context = std::dynamic_pointer_cast<GSContext<Mat,Vec>>(context_ptr_);
  gs_context->ApplyInverseTransportOperator();

  //=================================================== Assemble PETSc vector
  lbs_solver_.
  SetGSPETScVecFromPrimarySTLvector(groupset_,
                                    b_,
                                    lbs_solver_.PhiNewLocal(),
                                    WITH_DELAYED_PSI);

  //============================================= Compute RHS norm
  VecNorm(b_, NORM_2, &context_ptr_->rhs_norm_);

  //============================================= Compute precondition RHS norm
  PC pc;
  KSPGetPC(solver_, &pc);
  Vec temp_vec;
  VecDuplicate(b_, &temp_vec);
  PCApply(pc, b_, temp_vec);
  VecNorm(temp_vec, NORM_2, &context_ptr_->rhs_preconditioned_norm_);
}

/**Sets the initial guess for a gs solver. If the initial guess's norm
 * is large enough the initial guess will be used, otherwise it is assumed
 * zero.*/
template<> void GSLinearSolver<Mat, Vec, KSP>::SetInitialGuess()
{
  lbs_solver_.
  SetGSPETScVecFromPrimarySTLvector(groupset_, x_, lbs_solver_.PhiOldLocal());

  double init_guess_norm = 0.0;
  VecNorm(x_,NORM_2,&init_guess_norm);

  if (init_guess_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(solver_, PETSC_TRUE);
    if (log_info) chi::log.Log() << "Using phi_old as initial guess.";
  }
}

/**For this callback we simply restore the q_moments_local vector.*/
template<> void GSLinearSolver<Mat, Vec, KSP>::PostSolveCallback()
{
  lbs_solver_.QMomentsLocal() = saved_q_moments_local_;
}
}//namespace lbs

