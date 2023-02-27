#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

//###################################################################
/**Constructor.*/
chi_mesh::sweep_management::SweepBuffer::
  SweepBuffer(chi_mesh::sweep_management::AngleSet* ref_angleset,
              int sweep_eager_limit,
              const chi_objects::ChiMPICommunicatorSet& in_comm_set) :
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
bool chi_mesh::sweep_management::SweepBuffer::DoneSending() const
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

  const auto& spds =  angleset->GetSPDS();

  done_sending = true;
  for (size_t deplocI=0; deplocI<spds.location_successors.size(); deplocI++)
  {
    int num_mess = deplocI_message_count[deplocI];
    for (int m=0; m<num_mess; m++)
    {
      int message_sent = false;
      MPI_Test(&deplocI_message_request[deplocI][m],
               &message_sent, MPI_STATUS_IGNORE);
      if (not message_sent)
        done_sending = false;
    }
  }//for deplocI

  if (done_sending)
  {
    for (size_t deplocI=0; deplocI<spds.location_successors.size(); deplocI++)
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

//  for (auto & prelocI : prelocI_message_received)
//    for (auto && m : prelocI)
//      m = false;

  for (auto& message_flags : prelocI_message_received)
    message_flags.assign(message_flags.size(), false);

  for (auto& message_flags : delayed_prelocI_message_received)
    message_flags.assign(message_flags.size(), false);
}