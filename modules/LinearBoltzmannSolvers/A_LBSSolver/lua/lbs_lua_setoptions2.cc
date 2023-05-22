#include "lbs_lua_utils.h"

#include "ChiConsole/chi_console.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

namespace lbs::common_lua_utils
{

// ##################################################################
RegisterWrapperFunction(/*namespace_in_lua=*/lbs,
                        /*name_in_lua=*/SetOptions,
                        /*syntax_function=*/GetSyntax_SetOptions,
                        /*actual_function=*/SetOptions);

chi_objects::InputParameters GetSyntax_SetOptions()
{
  chi_objects::InputParameters params;

  // clang-format off
  params.SetGeneralDescription("\\defgroup lbs__SetOptions lbs.SetOptions \n"
                               "\\ingroup LBSLuaFunctions\n"
                               "Set options from a large list of parameters");

  params.AddRequiredParameter<size_t>(
    "arg0", "Handle to a <TT>lbs::LBSSolver</TT> object.");
  params.AddRequiredParameterBlock(
    "arg1", "Block of parameters $(lbs::OptionsBlock$)");

  // clang-format on

  return params;
}

chi_objects::ParameterBlock
SetOptions(const chi_objects::InputParameters& params)
{
  const std::string fname = __FUNCTION__;

  params.RequireParameter("arg0");
  params.RequireParameter("arg1");

  const size_t handle = params.GetParamValue<size_t>("arg0");
  auto& lbs_solver =
    chi::GetStackItem<lbs::LBSSolver>(chi::object_stack, handle, fname);

  auto options_params = LBSSolver::OptionsBlock();
  options_params.AssignParameters(params.GetParam("arg1"));

  lbs_solver.SetOptions(options_params);

  return chi_objects::ParameterBlock();
}

} // namespace lbs::common_lua_utils