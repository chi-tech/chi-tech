#ifndef CHITECH_RPPLOGICALVOLUME_H
#define CHITECH_RPPLOGICALVOLUME_H

#include "LogicalVolume.h"

namespace chi_mesh
{

// ###################################################################
/**Rectangular Parallel Piped (RPP) logical volume*/
class RPPLogicalVolume : public LogicalVolume
{
public:
  static chi_objects::InputParameters GetInputParameters();
  explicit RPPLogicalVolume(const chi_objects::InputParameters& params);

  bool Inside(const chi_mesh::Vector3& point) const override;

protected:
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;
};

} // namespace chi_mesh

#endif // CHITECH_RPPLOGICALVOLUME_H
