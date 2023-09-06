#ifndef CHITECH_FIELDFUNCTIONINTERFACE_H
#define CHITECH_FIELDFUNCTIONINTERFACE_H

#include "parameters/input_parameters.h"

namespace chi_physics
{

class FieldFunction;

/**A utility class to couple an object to a field function.*/
class FieldFunctionInterface
{
protected:
  static chi::InputParameters GetInputParameters();

  explicit FieldFunctionInterface(const chi::InputParameters& params);

  const FieldFunction* GetFieldFunction() const;

private:
  chi::ParameterBlock field_function_param_;
};

}

#endif // CHITECH_FIELDFUNCTIONINTERFACE_H
