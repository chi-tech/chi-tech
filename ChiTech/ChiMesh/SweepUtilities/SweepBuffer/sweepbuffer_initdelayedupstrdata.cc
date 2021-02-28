#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

//###################################################################
/**Initializes delayed upstream data. This method gets called
 * when a sweep scheduler is constructed.*/
void chi_mesh::sweep_management::SweepBuffer::
InitializeDelayedUpstreamData()
{
  auto  spds =  angleset->GetSPDS();
  auto fluds=  angleset->fluds;

  int num_grps   = angleset->GetNumGrps();
  int num_angles = angleset->angles.size();

  angleset->delayed_prelocI_outgoing_psi.clear();
  angleset->delayed_prelocI_outgoing_psi.resize(
    spds->delayed_location_dependencies.size());

  angleset->delayed_prelocI_outgoing_psi_old.clear();
  angleset->delayed_prelocI_outgoing_psi_old.resize(
    spds->delayed_location_dependencies.size());

  for (int prelocI=0;
       prelocI<spds->delayed_location_dependencies.size(); prelocI++)
  {
    int num_dofs = fluds->delayed_prelocI_face_dof_count[prelocI];

    u_ll_int buff_size = num_dofs*num_grps*num_angles;

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