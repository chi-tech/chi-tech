#include "../lbkes_k_eigenvalue_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

#include <petscsnes.h>

#define sc_int64_t static_cast<int64_t>

PetscErrorCode FormFunction(SNES, Vec, Vec, void *);

int lbs::DiscOrdKEigenvalueSolver::NonLinearKEigen()
{
  chi::log.Log()
    << "\n********** Solving k-eigenvalue problem with "
    << "Non-Linear iteration.\n";

  const auto& sdm = *discretization_;

  std::vector<int> groupset_ids;
  for (auto& groupset : groupsets_)
    groupset_ids.push_back(groupset.id_);

  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_globl_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info =
      groupset.angle_agg_.GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_globl_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_globl_dofs = num_globl_phi_dofs + num_globl_psi_dofs;

  chi::log.Log()
    << "Total number of angular unknowns: "
    << num_globl_psi_dofs;

  //============================================= Create the vectors
  Vec phi, r; /* solution, residual vectors */
  phi = chi_math::PETScUtils::CreateVector(sc_int64_t(num_local_dofs),
                                           sc_int64_t(num_globl_dofs));
  VecDuplicate(phi, &r);

  //============================================= Create SNES
  SNES        snes; /* nonlinear solver context */

  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetType(snes, SNESNEWTONLS);
  SNESSetApplicationContext(snes, &(*primary_ags_solver_->GetContext()));
  SNESSetFromOptions(snes);

  SNESSetTolerances(snes,
                    tolerance_,                        //absolute tolerance
                    1.0e-50,                           //relative tolerance
                    tolerance_,                        //solution tolerance
                    static_cast<int>(max_iterations_), //max iterations
                    -1);                               //max r-evals

  SNESSetMaxLinearSolveFailures(snes, 1000);

  //============================================= Set the residual function
  SNESSetFunction(snes, r, FormFunction, nullptr);

  //============================================= Setup Matrix-Free Jacobian
  Mat         J;    /* Jacobian matrix */

  MatCreateSNESMF(snes, &J);
  SNESSetJacobian(snes, J, J, MatMFFDComputeJacobian, nullptr);

  //============================================= Set linear solver and
  //                                              preconditioner
  KSP         ksp;  /* linear solver context */
  PC          pc;   /* preconditioner context */

  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);

  KSPSetTolerances(ksp,
                   1.0e-50,                 //relative tol
                   tolerance_,  //absolute tol
                   1.0e6,                   //divergence tol
//                   static_cast<int>(max_iterations_));     //max iterations
                   50);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCNONE);

  //============================================= Compute initial guess
  phi_old_local_.assign(num_local_phi_dofs, 1.0);

  DiscOrdSteadyStateSolver::SetMultiGSPETScVecFromPrimarySTLvector(
    groupset_ids, phi, PhiSTLOption::PHI_OLD);

  //============================================= Solve
  SNESSolve(snes, nullptr, phi);

  //============================================= Unpack solution
  SetPrimarySTLvectorFromGroupScopedPETScVec(
    groups_.front().id_, groups_.back().id_, phi, phi_old_local_);

  //============================================= Compute final k_eff
  k_eff_ = ComputeFissionProduction(phi_old_local_);

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log()
    << "        Final k-eigenvalue    :        "
    << std::setprecision(7) << k_eff_;
  chi::log.Log() << "\n";

  return 0;
}

PetscErrorCode FormFunction(SNES snes, Vec phi, Vec r, void*)
{
  lbs::AGSContext<Mat,Vec,KSP>* ags_context;
  SNESGetApplicationContext(snes, &ags_context);

  auto& lbs_solver = (lbs::DiscOrdKEigenvalueSolver&)ags_context->lbs_solver_;
  const auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& q_moments_local = lbs_solver.QMomentsLocal();

  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();
  auto primary_ags_solver = lbs_solver.GetPrimaryAGSSolver();

  const size_t num_local_dofs = q_moments_local.size();

  std::vector<int> groupset_ids;
  for (const auto& groupset : lbs_solver.Groupsets())
    groupset_ids.push_back(groupset.id_);

  //============================================= Disassemble phi vector
  lbs_solver.SetPrimarySTLvectorFromMultiGSPETScVecFrom(
    groupset_ids, phi, lbs::PhiSTLOption::PHI_OLD);

  //============================================= Compute 1/k F phi
  q_moments_local.assign(num_local_dofs, 0.0);
  for (auto& groupset : lbs_solver.Groupsets())
    active_set_source_function(groupset, q_moments_local, phi_old_local,
      lbs::APPLY_AGS_FISSION_SOURCES | lbs::APPLY_WGS_FISSION_SOURCES);

  const double k_eff = lbs_solver.ComputeFissionProduction(phi_old_local);
  chi_math::Scale(q_moments_local, 1.0/k_eff);

  //============================================= Now add MS phi
  for (auto& groupset : lbs_solver.Groupsets())
    active_set_source_function(groupset, q_moments_local, phi_old_local,
      lbs::APPLY_AGS_SCATTER_SOURCES | lbs::APPLY_WGS_SCATTER_SOURCES);

  //============================================= Sweep all the groupsets
  // After this phi_new = DLinv(MSD phi + 1/k FD phi)
  for (auto& groupset : lbs_solver.Groupsets())
  {
    auto& wgs_solver = lbs_solver.GetWGSSolvers()[groupset.id_];
    auto& raw_context = wgs_solver->GetContext();

    typedef lbs::WGSContext<Mat, Vec, KSP> LBSWGSContext;
    auto wgs_context = std::dynamic_pointer_cast<LBSWGSContext>(raw_context);

    wgs_context->ApplyInverseTransportOperator(lbs::NO_FLAGS_SET);
  }

  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(
    groupset_ids, r, lbs::PhiSTLOption::PHI_NEW);

  VecAXPY(r, -1.0, phi);

  return 0;
}