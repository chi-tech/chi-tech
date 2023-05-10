#include "lbs_DO_steady_state.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiObject/object_maker.h"

namespace lbs
{

RegisterChiObject(lbs,DiscOrdSteadyStateSolver);

// ###################################################################
/***/
chi_objects::InputParameters DiscOrdSteadyStateSolver::GetInputParameters()
{
  chi_objects::InputParameters params =
    DiscreteOrdinatesSolver::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "DiscOrdSteadyStateSolver");

  return params;
}

// ###################################################################
/**Static registration based constructor.*/
DiscOrdSteadyStateSolver::DiscOrdSteadyStateSolver(
  const chi_objects::InputParameters& params)
  : DiscreteOrdinatesSolver(params)
{
}

// ###################################################################
/**Execute the solver.*/
void DiscOrdSteadyStateSolver::Execute()
{
  primary_ags_solver_->Setup();
  primary_ags_solver_->Solve();

  if (options_.use_precursors) ComputePrecursors();

  UpdateFieldFunctions();

  chi::log.Log() << "LB solver " << TextName() << " execution completed\n";
}

} // namespace lbs
