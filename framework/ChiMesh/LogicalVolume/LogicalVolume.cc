#include "LogicalVolume.h"

namespace chi_mesh
{

chi_objects::InputParameters LogicalVolume::GetInputParameters()
{
  return ChiObject::GetInputParameters();
}

LogicalVolume::LogicalVolume(const chi_objects::InputParameters& params)
  : ChiObject(params)
{
}

} // namespace chi_mesh