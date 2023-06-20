#include "theta_scheme_time_intgr.h"

#include "ChiObject/object_maker.h"

#define scint static_cast<int>

namespace chi_math
{

RegisterChiObject(chi_math, ThetaSchemeTimeIntegration);

chi_objects::InputParameters ThetaSchemeTimeIntegration::GetInputParameters()
{
  chi_objects::InputParameters params = TimeIntegration::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_math__ThetaSchemeTimeIntegration "
  " chi_math.ThetaSchemeTimeIntegration\n"
  "\\ingroup DocTimeIntegrations");
  // clang-format on

  params.ChangeExistingParamToOptional("method",
                                       scint(SteppingMethod::THETA_SCHEME));

  params.AddRequiredParameter<double>("theta",
                                      "The theta parameter for a theta scheme");

  return params;
}

ThetaSchemeTimeIntegration::ThetaSchemeTimeIntegration(
  const chi_objects::InputParameters& params)
  : TimeIntegration(params), theta_(params.GetParamValue<double>("theta"))
{
}

double ThetaSchemeTimeIntegration::ThetaFactor() const { return theta_; }

} // namespace chi_math