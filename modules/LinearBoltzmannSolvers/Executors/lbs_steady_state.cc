#include "lbs_steady_state.h"

#include "ChiObject/object_maker.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

namespace lbs
{

RegisterChiObject(lbs, SteadyStateSolver);

chi_objects::InputParameters SteadyStateSolver::GetInputParameters()
{
  chi_objects::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription(
    "\\defgroup lbs__SteadyStateSolver lbs.SteadyStateSolver \n"
    "\\ingroup LBSExecutors\n"
    "Generalized implementation of a steady state solver. This solver calls"
    " the Across-Groupset (AGS) solver for the lbs-data block.");

  params.ChangeExistingParamToOptional("name", "SteadyStateSolver");

  params.AddRequiredParameter<size_t>("lbs_solver_handle",
                                      "Handle to an existing lbs solver");

  return params;
}

SteadyStateSolver::SteadyStateSolver(const chi_objects::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(chi::GetStackItem<LBSSolver>(
      chi::object_stack, params.GetParamValue<size_t>("lbs_solver_handle")))
{

}

void SteadyStateSolver::Initialize()
{
  lbs_solver_.Initialize();
}

void SteadyStateSolver::Execute()
{
  auto& ags_solver = *lbs_solver_.GetPrimaryAGSSolver();

  ags_solver.Setup();
  ags_solver.Solve();

  if (lbs_solver_.Options().use_precursors)
    lbs_solver_.ComputePrecursors();

  lbs_solver_.UpdateFieldFunctions();
}

} // namespace lbs