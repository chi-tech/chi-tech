#ifndef CHITECH_THETA_SCHEME_TIME_INTGR_H
#define CHITECH_THETA_SCHEME_TIME_INTGR_H

#include "time_integration.h"

namespace chi_math
{

class ThetaSchemeTimeIntegration : public TimeIntegration
{
private:
  const double theta_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit ThetaSchemeTimeIntegration(
    const chi_objects::InputParameters& params);

  double ThetaFactor() const;
};

} // namespace chi_math

#endif // CHITECH_THETA_SCHEME_TIME_INTGR_H
