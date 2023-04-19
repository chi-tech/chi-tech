#ifndef CHITECH_CHI_MATH_FUNCTION_H
#define CHITECH_CHI_MATH_FUNCTION_H

#include "ChiObject/chi_object.h"

namespace chi_math::functions
{

/**This is a base class for a complexity-limited function call. What
 * we want to support is the following.
 *   i)  Return of vector values (dim=1 to N).
 *   ii) Callable with
 *     a) Just a position
 *     b) Just the value of one or more variables
 *     c) With position AND and the value of one or more variables
 * and the ability.
 *
 * We can achieve this goal to some degree by defining the general function
\f[
F(x,y,z,[\phi_0, \dots, \phi_{N^a}]) \mapsto \mathbb{R}^b
\f]
where \f$ N^a \f$ is the dimension of the input values (this excludes the
position variables), and \f$ b \f$ is the number of output values.
*/
class FunctionXYZDimAToDimB : public ChiObject
{
private:
  const size_t input_dimension_;
  const size_t output_dimension_;
public:
  static chi_objects::InputParameters GetInputParameters();

  explicit FunctionXYZDimAToDimB(const chi_objects::InputParameters& params);

  /**Returns the input dimension.*/
  size_t InputDimension() const;
  /**Returns the output dimension.*/
  size_t OutputDimension() const;

  /**Evaluation at a point x.*/
  virtual std::vector<double> Evaluate(double x,
                                       double y,
                                       double z,
                                       const std::vector<double>& values) const;
};

} // namespace chi_math::functions

#endif // CHITECH_CHI_MATH_FUNCTION_H
