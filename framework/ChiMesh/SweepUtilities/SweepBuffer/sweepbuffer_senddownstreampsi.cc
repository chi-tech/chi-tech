#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include "ChiMPI/chi_mpi_commset.h"

//###################################################################
/**Sends downstream psi. This method gets called after a sweep chunk has
 * executed */
void chi_mesh::sweep_management::SweepBuffer::
SendDownstreamPsi(int angle_set_num)
{
  const auto& spds = angleset->GetSPDS();

  const size_t num_successors = spds.location_successors.size();
  for (size_t deplocI=0; deplocI<num_successors; deplocI++)
  {
    int locJ = spds.location_successors[deplocI];

    int num_mess = deplocI_message_count[deplocI];
    for (int m=0; m<num_mess; m++)
    {
      u_ll_int block_addr   = deplocI_message_blockpos[deplocI][m];
      u_ll_int message_size = deplocI_message_size[deplocI][m];

      const auto& outgoing_psi = angleset->deplocI_outgoing_psi[deplocI];

      MPI_Isend(&outgoing_psi[block_addr],
                static_cast<int>(message_size),
                MPI_DOUBLE,
                comm_set.MapIonJ(locJ,locJ),
                max_num_mess*angle_set_num + m, //tag
                comm_set.LocICommunicator(locJ),
                &deplocI_message_request[deplocI][m]);
    }//for message
  }//for deplocI
}