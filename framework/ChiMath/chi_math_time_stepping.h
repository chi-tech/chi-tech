#ifndef CHITECH_CHI_MATH_TIME_STEPPING_H
#define CHITECH_CHI_MATH_TIME_STEPPING_H

namespace chi_math
{

  enum class SteppingMethod
  {
    BACKWARD_EULER  = 0,
    CRANK_NICHOLSON = 1
  };

}//namespace chi_math

#endif //CHITECH_CHI_MATH_TIME_STEPPING_H
