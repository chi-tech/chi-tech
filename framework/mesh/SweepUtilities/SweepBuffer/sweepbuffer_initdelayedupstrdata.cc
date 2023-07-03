#include "sweepbuffer.h"

#include "mesh/SweepUtilities/AngleSet/angleset.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

//###################################################################
/**Initializes delayed upstream data. This method gets called
 * when a sweep scheduler is constructed.*/
void chi_mesh::sweep_management::SweepBuffer::
  InitializeDelayedUpstreamData()
{
  const auto& spds = angleset->GetSPDS();
  auto fluds=  angleset->fluds;

  const auto num_grps   = angleset->GetNumGrps();
  const auto num_angles = angleset->angles.size();

  const auto num_loc_deps = spds.GetDelayedLocationDependencies().size();

  angleset->delayed_prelocI_outgoing_psi.clear();
  angleset->delayed_prelocI_outgoing_psi.resize(
    num_loc_deps);

  angleset->delayed_prelocI_outgoing_psi_old.clear();
  angleset->delayed_prelocI_outgoing_psi_old.resize(
    num_loc_deps);

  for (int prelocI=0;prelocI<num_loc_deps; prelocI++)
  {
    const int num_nodes = fluds->delayed_prelocI_face_dof_count[prelocI];

    u_ll_int buff_size = num_nodes * num_grps * num_angles;

    angleset->delayed_prelocI_outgoing_psi[prelocI].resize(buff_size,0.0);
    angleset->delayed_prelocI_outgoing_psi_old[prelocI].resize(buff_size,0.0);
  }

  angleset->delayed_local_psi.resize(fluds->delayed_local_psi_stride*
                           fluds->delayed_local_psi_max_elements*
                           num_grps*num_angles,0.0);

  angleset->delayed_local_psi_old.resize(fluds->delayed_local_psi_stride*
                                         fluds->delayed_local_psi_max_elements*
                                         num_grps*num_angles,0.0);
}