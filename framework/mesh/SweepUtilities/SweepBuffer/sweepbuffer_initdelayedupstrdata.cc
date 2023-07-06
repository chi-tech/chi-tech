#include "sweepbuffer.h"

#include "mesh/SweepUtilities/AngleSet/angleset.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/AAH_FLUDS.h"

// ###################################################################
/**Initializes delayed upstream data. This method gets called
 * when a sweep scheduler is constructed.*/
void chi_mesh::sweep_management::SweepBuffer::InitializeDelayedUpstreamData()
{
  const auto& spds = angleset->GetSPDS();
  auto fluds = std::dynamic_pointer_cast<AAH_FLUDS>(angleset->fluds);

  const auto num_grps = angleset->GetNumGrps();
  const auto num_angles = angleset->angles.size();

  const auto num_loc_deps = spds.GetDelayedLocationDependencies().size();

  angleset->fluds->AllocateDelayedPrelocIOutgoingPsi(
    num_grps, num_angles, num_loc_deps);

  angleset->fluds->AllocateDelayedLocalPsi(num_grps, num_angles);
}