#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

//###################################################################
/**Sends downstream psi. This method gets called after a sweep chunk has
 * executed */
void chi_mesh::sweep_management::SweepBuffer::
SendDownstreamPsi(int angle_set_num)
{
  chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();

  for (size_t deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
  {
    int locJ = spds->location_successors[deplocI];

    int num_mess = deplocI_message_count[deplocI];
    for (int m=0; m<num_mess; m++)
    {
      u_ll_int block_addr   = deplocI_message_blockpos[deplocI][m];
      u_ll_int message_size = deplocI_message_size[deplocI][m];

      MPI_Isend(&angleset->deplocI_outgoing_psi[deplocI].data()[block_addr],
                message_size,
                MPI_DOUBLE,
                comm_set->MapIonJ(locJ,locJ),
                max_num_mess*angle_set_num + m, //tag
                comm_set->communicators[locJ],
                &deplocI_message_request[deplocI][m]);
    }//for message
  }//for deplocI
}