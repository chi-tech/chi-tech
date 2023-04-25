#include "xxpoweriteration_keigen.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include <iomanip>

namespace lbs
{

RegisterChiObject(lbs, XXPowerIterationKEigen);

chi_objects::InputParameters XXPowerIterationKEigen::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

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

  return params;
}

XXPowerIterationKEigen::XXPowerIterationKEigen(
  const chi_objects::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(chi::GetStackItem<LBSSolver>(
      chi::object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    max_iters_(params.GetParamValue<size_t>("max_iters")),
    k_tolerance_(params.GetParamValue<double>("k_tol")),
    reset_solution_(params.GetParamValue<bool>("reset_solution")),

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

    wgs_context->lhs_src_scope_ = APPLY_WGS_SCATTER_SOURCES;
    wgs_context->rhs_src_scope_ =
      APPLY_AGS_SCATTER_SOURCES | APPLY_FIXED_SOURCES;
  }

  primary_ags_solver_.SetVerbosity(
    lbs_solver_.Options().verbose_ags_iterations);

  front_wgs_solver_ = lbs_solver_.GetWGSSolvers().at(front_gs_.id_);
  front_wgs_context_ =
    std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(
      front_wgs_solver_->GetContext());

  ChiLogicalErrorIf(not front_wgs_context_, fname + ": Casting failure");

  if (reset_solution_)
    lbs_solver_.SetPhiVectorScalarValues(phi_old_local_, 1.0);
}

// ##################################################################
/**Empty initialization function.*/
void XXPowerIterationKEigen::Initialize() {}

// ##################################################################
/**Executes the solver.*/
void XXPowerIterationKEigen::Execute()
{
  using namespace chi_math;

  double F_prev = 1.0;
  k_eff_ = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iters_)
  {
    //================================= Set the fission source
    SetLBSFissionSource(phi_old_local_, /*additive=*/false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    //================================= This solves the inners for transport
    primary_ags_solver_.Setup();
    primary_ags_solver_.Solve();

    //================================= Recompute k-eigenvalue
    double F_new = lbs_solver_.ComputeFissionProduction(phi_new_local_);
    k_eff_ = F_new / F_prev * k_eff_;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    //================================= Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(k_tolerance_, 1.0e-12)) converged = true;

    //================================= Print iteration summary
    if (lbs_solver_.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << chi::program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff "
                  << std::setw(11) << std::setprecision(7) << k_eff_
                  << "  k_eff change " << std::setw(12) << k_eff_change
                  << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  } // for k iterations

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log() << "        Final k-eigenvalue    :        "
                 << std::setprecision(7) << k_eff_;
  chi::log.Log() << "        Final change          :        "
                 << std::setprecision(6) << k_eff_change << " (num_TrOps:"
                 << front_wgs_context_->counter_applications_of_inv_op_ << ")"
                 << "\n";
  chi::log.Log() << "\n";

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    chi_math::Scale(lbs_solver_.PrecursorsNewLocal(), 1.0 / k_eff_);
  }

  lbs_solver_.UpdateFieldFunctions();

  chi::log.Log()
    << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

// ##################################################################
/**Combines function calls to set fission source.*/
void XXPowerIterationKEigen::SetLBSFissionSource(const VecDbl& input,
                                                 const bool additive)
{
  if (not additive) chi_math::Set(q_moments_local_, 0.0);
  active_set_source_function_(front_gs_,
                              q_moments_local_,
                              input,
                              APPLY_AGS_FISSION_SOURCES |
                                APPLY_WGS_FISSION_SOURCES);
}

// ##################################################################
/**Combines function calls to set scattering source source.*/
void XXPowerIterationKEigen::SetLBSScatterSource(
  const VecDbl& input, const bool additive, const bool suppress_wgs /*=false*/)
{
  if (not additive) chi_math::Set(q_moments_local_, 0.0);
  active_set_source_function_(
    front_gs_,
    q_moments_local_,
    input,
    APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
      (suppress_wgs ? SUPPRESS_WG_SCATTER : NO_FLAGS_SET));
}

} // namespace lbs