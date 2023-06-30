#include "sweepbuffer.h"

#include "mesh/SweepUtilities/AngleSet/angleset.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

#include "chi_runtime.h"
#include "console/chi_console.h"
#include "chi_log.h"

//###################################################################
/** This is the final level of initialization before a sweep-chunk executes.
 * Once all upstream dependencies are met and if the sweep scheduler places
 * this angleset as "ready-to-execute", then the angle-set will call this
 * method. It is also fairly important in terms of memory to only allocate
 * these chunks of memory when actually ready to use them since they form the
 * majority of memory usage.*/
void chi_mesh::sweep_management::SweepBuffer::
  InitializeLocalAndDownstreamBuffers()
{
  if (!data_initialized)
  {
    const auto& spds = angleset->GetSPDS();
    auto fluds=  angleset->fluds;

    const auto num_grps   = angleset->GetNumGrps();
    const auto num_angles = angleset->angles.size();

    //============================ Resize FLUDS local outgoing Data
    angleset->local_psi.resize(fluds->num_face_categories);
    // fc = face category
    for (size_t fc = 0; fc<fluds->num_face_categories; fc++)
    {
      angleset->local_psi[fc].resize(fluds->local_psi_stride[fc]*
                                     fluds->local_psi_max_elements[fc]*
                                     num_grps*num_angles,0.0);
    }

    //============================ Resize FLUDS non-local outgoing Data
    angleset->deplocI_outgoing_psi.resize(
      spds.location_successors.size(),std::vector<double>());
    for (size_t deplocI=0; deplocI<spds.location_successors.size(); deplocI++)
    {
      angleset->deplocI_outgoing_psi[deplocI].resize(
        fluds->deplocI_face_dof_count[deplocI]*num_grps*num_angles,0.0);
    }

    //================================================ Make a memory query
    double memory_mb = chi::ChiConsole::GetMemoryUsageInMB();

    std::shared_ptr<chi::ChiLog::EventInfo> memory_event_info =
      std::make_shared<chi::ChiLog::EventInfo>(memory_mb);

    Chi::log.LogEvent(chi::ChiLog::StdTags::MAX_MEMORY_USAGE,
                      chi::ChiLog::EventType::SINGLE_OCCURRENCE,
                     memory_event_info);

    data_initialized = true;
  }
}