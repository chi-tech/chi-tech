#include "chi_angleset.h"
#include "chi_SPDS.h"
#include "chi_sweepbuffer.h"

#include <chi_mpi.h>

extern ChiMPI chi_mpi;

//###################################################################
/**Initializes delayed upstream data.*/
void chi_mesh::SweepManagement::AngleSet::
InitializeDelayedUpstreamData()
{
  int num_angles = angles.size();

  delayed_prelocI_outgoing_psi.resize(
      spds->delayed_location_dependencies.size(),
      std::vector<double>());
  for (int prelocI=0;
       prelocI<spds->delayed_location_dependencies.size(); prelocI++)
  {
    int num_dofs = fluds->delayed_prelocI_face_dof_count[prelocI];
    int num_grps   = GetNumGrps();

    u_ll_int buff_size = num_dofs*num_grps*num_angles;

    delayed_prelocI_outgoing_psi[prelocI].resize(buff_size,0.0);
  }

  delayed_local_psi.resize(fluds->delayed_local_psi_stride*
                           fluds->delayed_local_psi_max_elements*
                           num_grps*num_angles,0.0);
}

//###################################################################
/**This function advances the work stages of an angleset.*/
bool chi_mesh::SweepManagement::AngleSet::
AngleSetAdvance(chi_mesh::SweepManagement::SweepChunk *sweep_chunk,
                int angle_set_num)
{
  //================================================== Prevent reexecution
  if (executed)
    return FLAG_FINISHED;

  //================================================== Check all predecessor
  //                                                   locations sent data
  if (!sweep_buffer.CheckUpstreamPsiAvailable(angle_set_num))
    return FLAG_NOT_FINISHED;

  //================================================== Resize the angle set's
  //                                                   FLUDS at workstage start
  sweep_buffer.InitializeBuffers();

  //================================================== Execute chunk
  sweep_chunk->Sweep(this);

  //================================================== Send outgoing psi
  sweep_buffer.SendDownstreamPsi(angle_set_num);

  //================================================== Release memory
  sweep_buffer.ClearReceiveBuffers();

  executed = true;
  return FLAG_FINISHED;
}