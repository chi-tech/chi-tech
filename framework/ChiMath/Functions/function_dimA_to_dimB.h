#ifndef CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H
#define CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H

#include "ChiObject/chi_object.h"

namespace chi_math
{

class FunctionDimAToDimB : public ChiObject
{
private:
  const size_t input_dimension_;
  const size_t output_dimension_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit FunctionDimAToDimB(const chi_objects::InputParameters& params);

  virtual std::vector<double>
  Evaluate(const std::vector<double>& vals) const = 0;
};

} // namespace chi_math

#endif // CHITECH_CHI_MATH_FUNCTION_DIMA_TO_DIMB_H
