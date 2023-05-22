#ifndef CHITECH_FIELD_OPERATION_H
#define CHITECH_FIELD_OPERATION_H

#include "ChiObject/chi_object.h"

namespace chi_physics::field_operations
{

/**The base field operation class.*/
class FieldOperation : public ChiObject
{
public:
  static chi_objects::InputParameters GetInputParameters();

  explicit FieldOperation(const chi_objects::InputParameters& params);

  virtual void Execute() = 0;

  virtual ~FieldOperation() = default;
};

}

#endif // CHITECH_FIELD_OPERATION_H
