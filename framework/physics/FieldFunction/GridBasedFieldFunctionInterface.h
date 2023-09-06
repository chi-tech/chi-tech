#ifndef CHITECH_GRIDBASEDFIELDFUNCTIONINTERFACE_H
#define CHITECH_GRIDBASEDFIELDFUNCTIONINTERFACE_H

#include "FieldFunctionInterface.h"

namespace chi_physics
{

class FieldFunctionGridBased;

class GridBasedFieldFunctionInterface : public FieldFunctionInterface
{
public:
  static chi::InputParameters GetInputParameters();

  explicit GridBasedFieldFunctionInterface(const chi::InputParameters& params);

  const FieldFunctionGridBased* GetGridBasedFieldFunction() const;
};

}

#endif // CHITECH_GRIDBASEDFIELDFUNCTIONINTERFACE_H
