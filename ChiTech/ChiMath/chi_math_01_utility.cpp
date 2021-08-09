#include "chi_math.h"

//###################################################################
/**Computes the factorial of an integer.*/
double chi_math::Factorial(const int x)
{
  double factorial_value = 1.0;
  for (int i=2; i<=x; ++i)
    factorial_value *= i;

  return factorial_value;
}



