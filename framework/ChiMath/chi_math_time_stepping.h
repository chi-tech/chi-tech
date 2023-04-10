#ifndef CHITECH_CHI_MATH_TIME_STEPPING_H
#define CHITECH_CHI_MATH_TIME_STEPPING_H

#include <string>

namespace chi_math
{

enum class SteppingMethod
{
  EXPLICIT_EULER = 0,
  IMPLICIT_EULER = 1,
  CRANK_NICOLSON = 2
};

std::string SteppingMethodStringName(SteppingMethod);
SteppingMethod SteppingMethodFromString(const std::string& name);

} // namespace chi_math

#endif // CHITECH_CHI_MATH_TIME_STEPPING_H
