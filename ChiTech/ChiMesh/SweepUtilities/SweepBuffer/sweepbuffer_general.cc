#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

//###################################################################
/**Constructor.*/
chi_mesh::sweep_management::SweepBuffer::
SweepBuffer(chi_mesh::sweep_management::AngleSet* ref_angleset,
            int sweep_eager_limit,
            chi_objects::ChiMPICommunicatorSet* in_comm_set):
  angleset(ref_angleset),
  comm_set(in_comm_set)
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;
  EAGER_LIMIT = sweep_eager_limit;

  max_num_mess = 0;
}

//###################################################################
/**Returns the private flag done_sending.*/
bool chi_mesh::sweep_management::SweepBuffer::DoneSending()
{
  return done_sending;
}

//###################################################################
/**Receive all upstream Psi. This method is called from within
 * an advancement of an angleset, right after execution.*/
void chi_mesh::sweep_management::SweepBuffer::
ClearLocalAndReceiveBuffers()
{
  auto empty_vector = std::vector<std::vector<double>>(0);
  angleset->local_psi.swap(empty_vector);

  empty_vector = std::vector<std::vector<double>>(0);
  angleset->prelocI_outgoing_psi.swap(empty_vector);
}

//###################################################################
/**Sends downstream psi.*/
void chi_mesh::sweep_management::SweepBuffer::
ClearDownstreamBuffers()
{
  if (done_sending) return;

  auto spds =  angleset->GetSPDS();

  done_sending = true;
  for (size_t deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
  {
    int num_mess = deplocI_message_count[deplocI];
    for (int m=0; m<num_mess; m++)
    {
      int  send_request_status = 1;
      MPI_Test(&deplocI_message_request[deplocI][m],
               &send_request_status,MPI_STATUS_IGNORE);
      if (send_request_status == 0) done_sending = false;
    }

  }

  if (done_sending)
  {
    for (size_t deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
    {
      angleset->deplocI_outgoing_psi[deplocI].clear();
      angleset->deplocI_outgoing_psi[deplocI].shrink_to_fit();
    }
  }
}

//###################################################################
/**Reset flags in preperation for another sweep.*/
void chi_mesh::sweep_management::SweepBuffer::Reset()
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;

  for (int prelocI=0; prelocI<prelocI_message_available.size(); prelocI++)
    for (int m=0; m<prelocI_message_available[prelocI].size(); m++)
      prelocI_message_available[prelocI][m] = false;

  for (int prelocI=0; prelocI<delayed_prelocI_message_available.size(); prelocI++)
    for (int m=0; m<delayed_prelocI_message_available[prelocI].size(); m++)
      delayed_prelocI_message_available[prelocI][m] = false;

}