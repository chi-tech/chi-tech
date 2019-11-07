#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

#include <ChiConsole/chi_console.h>

extern ChiConsole chi_console;

extern double chi_global_timings[20];

//###################################################################
/** This is the final level of initialization before a sweep-chunk executes.
 * Once all upstream dependencies are met and if the sweep scheduler places
 * this angleset as "ready-to-execute", then the angle-set will call this
 * method. It is also fairly important in terms of memory to only allocate
 * these chunks of memory since they form the majority of memory requirements.*/
void chi_mesh::sweep_management::SweepBuffer::
InitializeLocalAndDownstreamBuffers()
{
  if (!data_initialized)
  {
    chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();
    chi_mesh::sweep_management::FLUDS* fluds=  angleset->fluds;

    int num_grps   = angleset->GetNumGrps();
    int num_angles = angleset->angles.size();

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
      spds->location_successors.size(),std::vector<double>());
    for (size_t deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
    {
      angleset->deplocI_outgoing_psi[deplocI].resize(
        fluds->deplocI_face_dof_count[deplocI]*num_grps*num_angles,0.0);
    }

    //================================================ Make a memory query
    double memory_mb = chi_console.GetMemoryUsageInMB();

    if (memory_mb>chi_global_timings[9])
      chi_global_timings[9] = memory_mb;

    chi_global_timings[10] += memory_mb;
    chi_global_timings[11] += 1.0;

    data_initialized = true;
  }

  //================================================== Copy delayed Psi to Psi_old
  angleset->delayed_local_psi_old.clear();
  std::copy(angleset->delayed_local_psi.begin(),
            angleset->delayed_local_psi.end(),
            std::back_inserter(angleset->delayed_local_psi_old));
}