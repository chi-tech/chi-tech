#include "crank_nicolson_time_intgr.h"

#include "ChiObject/object_maker.h"

#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, CrankNicolsonTimeIntegration);

chi::InputParameters CrankNicolsonTimeIntegration::GetInputParameters()
{
  chi::InputParameters params =
    ThetaSchemeTimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_math__CrankNicolsonTimeIntegration "
  " chi_math.CrankNicolsonTimeIntegration\n"
  "\\ingroup DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method",
                                       scint(SteppingMethod::CRANK_NICOLSON));
  params.ChangeExistingParamToOptional("theta", 0.5);

  return params;
}

CrankNicolsonTimeIntegration::CrankNicolsonTimeIntegration(
  const chi::InputParameters& params)
  : ThetaSchemeTimeIntegration(params)
{
}

} // namespace chi_math