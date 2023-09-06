#ifndef CHITECH_LOGICALVOLUMEINTERFACE_H
#define CHITECH_LOGICALVOLUMEINTERFACE_H

#include "parameters/input_parameters.h"

namespace chi_mesh
{

class LogicalVolume;

/**Interface class to add a dependency on a logical volume.*/
class LogicalVolumeInterface
{
protected:
  static chi::InputParameters GetInputParameters();

  explicit LogicalVolumeInterface(const chi::InputParameters& params);

  const LogicalVolume* GetLogicalVolume() const;

private:
  const std::shared_ptr<const LogicalVolume> logical_volume_;
};

}

#endif // CHITECH_LOGICALVOLUMEINTERFACE_H
