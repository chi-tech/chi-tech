#include "chi_math_time_stepping.h"

#include <stdexcept>

namespace chi_math
{

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wswitch-enum"
/**Returns the string name of a time stepping method.*/
std::string SteppingMethodStringName(SteppingMethod method)
{
  switch (method)
  {
    case SteppingMethod::EXPLICIT_EULER: return "explicit_euler";
    case SteppingMethod::IMPLICIT_EULER: return "implicit_euler";
    case SteppingMethod::CRANK_NICOLSON: return "crank_nicholson";
    default:
      throw std::logic_error(__PRETTY_FUNCTION__);
  }
}
#pragma GCC diagnostic pop

}