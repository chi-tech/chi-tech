#ifndef CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H
#define CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H

#include "functionxyz_dimA_to_dimB.h"

namespace chi_math::functions
{
class PiecewiseLinear1D : public FunctionXYZDimAToDimB
{
private:
  const std::vector<double> x_values_;
  const std::vector<double> y_values_;
  const size_t num_vals_;
  std::vector<double> delta_x_values_;

public:
  static chi_objects::InputParameters GetInputParameters();

  explicit PiecewiseLinear1D(const chi_objects::InputParameters& params);

  /**Function evaluation at a single value.*/
  std::vector<double>
  Evaluate(double x,
           double y,
           double z,
           const std::vector<double>& values) const override;
};
} // namespace chi_math::functions

#endif // CHITECH_CHI_MATH_FUNCTIONS_PIECEWISELINEAR1D_H
