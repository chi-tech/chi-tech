#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/SweepBuffer/sweepbuffer.h"

#include <chi_mpi.h>

extern ChiMPI chi_mpi;



//###################################################################
/**Initializes delayed upstream data.*/
void chi_mesh::sweep_management::AngleSet::
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
chi_mesh::sweep_management::AngleSetStatus
chi_mesh::sweep_management::AngleSet::
AngleSetAdvance(chi_mesh::sweep_management::SweepChunk *sweep_chunk,
                int angle_set_num,
                chi_mesh::sweep_management::ExecutionPermission permission)
{
  typedef AngleSetStatus Status;

  if (executed)
  {
    if (!sweep_buffer.done_sending)
      sweep_buffer.ClearDownstreamBuffers();
    return AngleSetStatus::FINISHED;
  }

  Status status = sweep_buffer.ReceiveUpstreamPsi(angle_set_num);

  if      (status == Status::RECEIVING) return status;
  else if (status == Status::READY_TO_EXECUTE and
           permission == ExecutionPermission::EXECUTE)
  {
    sweep_buffer.InitializeBuffers();

    sweep_chunk->Sweep(this); //Execute chunk

    //Send outgoing psi and clear local and receive buffers
    sweep_buffer.SendDownstreamPsi(angle_set_num);
    sweep_buffer.ClearLocalAndReceiveBuffers();

    executed = true;
    return AngleSetStatus::FINISHED;
  }
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}