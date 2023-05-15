#include "pi_keigen.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

namespace lbs
{

RegisterChiObject(lbs, XXPowerIterationKEigen);

chi_objects::InputParameters XXPowerIterationKEigen::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription(
    "\\defgroup lbs__XXPowerIterationKEigen lbs.XXPowerIterationKEigen \n"
    "\\ingroup lbs__LBSSolver\n"
    "Generalized implementation of a k-Eigenvalue solver using Power "
    "Iteration.");

  params.ChangeExistingParamToOptional("name", "XXPowerIterationKEigen");

  params.AddRequiredParameter<size_t>("lbs_solver_handle",
                                      "Handle to an existing lbs solver");
  params.AddOptionalParameter(
    "max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter(
    "k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution",
    true,
    "Flag, if set to true will initialize the phi-solution to all 1's before "
    "executing");

  params.AddOptionalParameter(
    "reinit_phi_1", true, "If true, reinitializes scalar phi fluxes to 1");

  return params;
}

XXPowerIterationKEigen::XXPowerIterationKEigen(
  const chi_objects::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(chi::GetStackItem<LBSSolver>(
      chi::object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    max_iters_(params.GetParamValue<size_t>("max_iters")),
    k_tolerance_(params.GetParamValue<double>("k_tol")),
    reinit_phi_1_(params.GetParamValue<bool>("reinit_phi_1")),

    q_moments_local_(lbs_solver_.QMomentsLocal()),
    phi_old_local_(lbs_solver_.PhiOldLocal()),
    phi_new_local_(lbs_solver_.PhiNewLocal()),
    groupsets_(lbs_solver_.Groupsets()),
    primary_ags_solver_(*lbs_solver_.GetPrimaryAGSSolver()),
    active_set_source_function_(lbs_solver_.GetActiveSetSourceFunction()),
    front_gs_(groupsets_.front())
{
  const std::string fname = __PRETTY_FUNCTION__;

  for (auto& wgs_solver : lbs_solver_.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context =
      std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(context);

    ChiLogicalErrorIf(not wgs_context, fname + ": Cast failed");

    wgs_context->lhs_src_scope_ =
      wgs_context->lhs_src_scope_ & (~APPLY_WGS_FISSION_SOURCES); // lhs_scope
    wgs_context->rhs_src_scope_ =
      wgs_context->rhs_src_scope_ & (~APPLY_AGS_FISSION_SOURCES); // rhs_scope
  }

  primary_ags_solver_.SetVerbosity(
    lbs_solver_.Options().verbose_ags_iterations);

  front_wgs_solver_ = lbs_solver_.GetWGSSolvers().at(front_gs_.id_);
  front_wgs_context_ =
    std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(
      front_wgs_solver_->GetContext());

  ChiLogicalErrorIf(not front_wgs_context_, fname + ": Casting failure");

  if (reinit_phi_1_) lbs_solver_.SetPhiVectorScalarValues(phi_old_local_, 1.0);
}

} // namespace lbs