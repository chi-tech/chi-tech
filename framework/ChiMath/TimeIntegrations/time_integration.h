#ifndef CHITECH_TIME_INTEGRATION_H
#define CHITECH_TIME_INTEGRATION_H

#include "ChiObject/chi_object.h"
#include "ChiMath/chi_math_time_stepping.h"

namespace chi_math
{

class TimeIntegration : public ChiObject
{
private:
  SteppingMethod method_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit TimeIntegration(const chi_objects::InputParameters& params);

  SteppingMethod Method() const;

  virtual ~TimeIntegration() = default;
};

} // namespace chi_math

#endif // CHITECH_TIME_INTEGRATION_H
