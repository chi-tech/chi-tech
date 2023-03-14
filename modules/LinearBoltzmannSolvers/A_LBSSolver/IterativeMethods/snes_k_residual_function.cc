#include "A_LBSSolver/lbs_solver.h"
#include "A_LBSSolver/IterativeMethods/ags_context.h"
#include "A_LBSSolver/IterativeMethods/wgs_context.h"
#include "A_LBSSolver/Preconditioning/lbs_shell_operations.h"

#include "snes_k_residual_func_context.h"

#include <petscsnes.h>

namespace lbs
{

/**This function evaluates the flux moments based k-eigenvalue transport
 * residual of the form
\f$ r(\phi) = DL^{-1} (\frac{1}{k} F\phi + MS \phi) - \phi \f$.*/
PetscErrorCode SNESKResidualFunction(SNES snes, Vec phi, Vec r, void* ctx)
{
  const std::string fname  = "lbs::SNESKResidualFunction";
  auto& function_context = *((KResidualFunctionContext*)ctx);

  lbs::AGSContext<Mat,Vec,KSP>* ags_context;
  SNESGetApplicationContext(snes, &ags_context);

  auto& lbs_solver = ags_context->lbs_solver_;
  const auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& q_moments_local = lbs_solver.QMomentsLocal();

  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();

  std::vector<int> groupset_ids;
  for (const auto& groupset : lbs_solver.Groupsets())
    groupset_ids.push_back(groupset.id_);

  //============================================= Disassemble phi vector
  lbs_solver.SetPrimarySTLvectorFromMultiGSPETScVecFrom(
    groupset_ids, phi, lbs::PhiSTLOption::PHI_OLD);

  //============================================= Compute 1/k F phi
  chi_math::Set(q_moments_local, 0.0);
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

  //============================================= Reassemble PETSc vector
  //We use r as a proxy for delta-phi here since
  // we are anycase going to subtract phi from it.
  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(
    groupset_ids, r, lbs::PhiSTLOption::PHI_NEW);

  VecAXPY(r, -1.0, phi);

  for (auto& groupset : lbs_solver.Groupsets())
  {
    if ((groupset.apply_wgdsa_ or groupset.apply_tgdsa_) and
        lbs_solver.Groupsets().size() > 1)
      throw std::logic_error(fname + ": Preconditioning currently only supports"
                                     "single groupset simulations.");
    auto& wgs_solver = lbs_solver.GetWGSSolvers()[groupset.id_];
    auto& raw_context = wgs_solver->GetContext();
    auto wgs_context =
      std::dynamic_pointer_cast<lbs::WGSContext<Mat,Vec,KSP>>(raw_context);
    lbs::WGDSA_TGDSA_PreConditionerMult2(*wgs_context, r, r);
  }

  //============================================= Assign k to the context
  //                                              so monitors can work
  function_context.k_eff = k_eff;

  return 0;
}

}//namespace lbs