#include "nl_keigen.h"

#include "ChiObject/object_maker.h"
#include "chi_log.h"

namespace lbs
{

RegisterChiObject(lbs, XXNonLinearKEigen);

chi_objects::InputParameters XXNonLinearKEigen::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription(
    "\\defgroup lbs__XXNonLinearKEigen lbs.XXNonLinearKEigen \n"
    "\\ingroup LBSExecutors\n"
    "Generalized implementation of a non-linear k-Eigenvalue solver");

  params.ChangeExistingParamToOptional("name", "XXPowerIterationKEigen");

  params.AddRequiredParameter<size_t>("lbs_solver_handle",
                                      "Handle to an existing lbs solver");

  params.AddOptionalParameter(
    "nl_abs_tol", 1.0e-8, "Non-linear absolute tolerance");
  params.AddOptionalParameter(
    "nl_rel_tol", 1.0e-8, "Non-linear relative tolerance");
  params.AddOptionalParameter(
    "nl_sol_tol", 1.0e-50, "Non-linear solution tolerance");
  params.AddOptionalParameter(
    "nl_max_its", 50, "Non-linear maximum iterations");

  params.AddOptionalParameter("l_abs_tol", 1.0e-8, "Linear absolute tolerance");
  params.AddOptionalParameter("l_rel_tol", 1.0e-8, "Linear relative tolerance");
  params.AddOptionalParameter(
    "l_div_tol", 1.0e6, "Linear divergence tolerance");
  params.AddOptionalParameter("l_max_its", 50, "Linear maximum iterations");
  params.AddOptionalParameter(
    "l_gmres_restart_intvl", 30, "GMRes restart interval");
  params.AddOptionalParameter(
    "l_gmres_breakdown_tol", 1.0e6, "GMRes breakdown tolerance");

  params.AddOptionalParameter(
    "reinit_phi_1", true, "If true, reinitializes scalar phi fluxes to 1");

  return params;
}

XXNonLinearKEigen::XXNonLinearKEigen(const chi_objects::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(chi::GetStackItem<LBSSolver>(
      chi::object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    nl_context_(std::make_shared<NLKEigenAGSContext<Vec, SNES>>(lbs_solver_)),
    nl_solver_(SNESNEWTONLS, nl_context_),
    reinit_phi_1_(params.GetParamValue<bool>("reinit_phi_1"))
{
  auto& tolerances = nl_solver_.ToleranceOptions();

  tolerances.nl_absolute_tol = params.GetParamValue<double>("nl_abs_tol");
  tolerances.nl_relative_tol = params.GetParamValue<double>("nl_rel_tol");
  tolerances.nl_solution_tol = params.GetParamValue<double>("nl_sol_tol");
  tolerances.nl_max_iterations = params.GetParamValue<int>("nl_max_its");

  tolerances.l_relative_tol = params.GetParamValue<double>("l_rel_tol");
  tolerances.l_absolute_tol = params.GetParamValue<double>("l_abs_tol");
  tolerances.l_divergence_tol = params.GetParamValue<double>("l_div_tol");
  tolerances.l_max_iterations = params.GetParamValue<int>("l_max_its");
  tolerances.l_gmres_restart_interval =
    params.GetParamValue<int>("l_gmres_restart_intvl");
  tolerances.l_gmres_breakdown_tol =
    params.GetParamValue<double>("l_gmres_breakdown_tol");
}

void XXNonLinearKEigen::Initialize()
{
  lbs_solver_.Initialize();
}

void XXNonLinearKEigen::Execute()
{
  if (reinit_phi_1_)
    lbs_solver_.SetPhiVectorScalarValues(lbs_solver_.PhiOldLocal(), 1.0);

  nl_solver_.Setup();
  nl_solver_.Solve();

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    chi_math::Scale(lbs_solver_.PrecursorsNewLocal(),
                    1.0 / nl_context_->kresid_func_context_.k_eff);
  }

  lbs_solver_.UpdateFieldFunctions();

  chi::log.Log()
    << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

} // namespace lbs