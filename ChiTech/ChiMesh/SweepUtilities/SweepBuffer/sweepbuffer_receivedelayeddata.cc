#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/** Receives delayed data from successor locations. */
bool chi_mesh::sweep_management::SweepBuffer::
  ReceiveDelayedData(int angle_set_num)
{
  auto spds =  angleset->GetSPDS();

  const size_t num_delayed_loc_deps = spds->delayed_location_dependencies.size();

  //======================================== Receive delayed data
  bool all_messages_received = true;
  for (size_t prelocI=0; prelocI<num_delayed_loc_deps; prelocI++)
  {
    int locJ = spds->delayed_location_dependencies[prelocI];


    int num_mess = delayed_prelocI_message_count[prelocI];
    for (int m=0; m<num_mess; m++)
    {
      if (not delayed_prelocI_message_received[prelocI][m])
      {
        int message_available = 0;
        MPI_Iprobe(comm_set->MapIonJ(locJ,chi::mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi::mpi.location_id],
                   &message_available, MPI_STATUS_IGNORE);

        if (not message_available)
        {
          all_messages_received = false;
          continue;
        }

        //============================ Receive upstream data
        auto& upstream_psi = angleset->delayed_prelocI_outgoing_psi[prelocI];

        u_ll_int block_addr   = delayed_prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = delayed_prelocI_message_size[prelocI][m];

        int error_code =
          MPI_Recv(&upstream_psi[block_addr],
                   static_cast<int>(message_size),
                   MPI_DOUBLE,
                   comm_set->MapIonJ(locJ,chi::mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi::mpi.location_id],
                   MPI_STATUS_IGNORE);

        delayed_prelocI_message_received[prelocI][m] = true;

        if (error_code != MPI_SUCCESS)
        {
          std::stringstream err_stream;
          err_stream << "################# Delayed receive error."
                     << " message size=" << message_size
                     << " as_num=" << angle_set_num
                     << " num_mess=" << num_mess
                     << " m=" << m
                     << " error="
                     << " size=" << "\n";
          char error_string[BUFSIZ];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          MPI_Error_string(error_code, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          chi::log.LogAllWarning() << err_stream.str();
        }
      }//if not message already received
    }//for message
  }//for delayed predecessor

  if (not all_messages_received)
    return false;

  return true;
}