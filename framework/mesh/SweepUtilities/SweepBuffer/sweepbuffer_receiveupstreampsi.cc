#include "sweepbuffer.h"

#include "mesh/SweepUtilities/AngleSet/angleset.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

#include "mpi/chi_mpi_commset.h"

#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Check if all upstream dependencies have been met and receives
 * it as it becomes available.*/
chi_mesh::sweep_management::AngleSetStatus
chi_mesh::sweep_management::SweepBuffer::ReceiveUpstreamPsi(int angle_set_num)
{
  const auto& spds = angleset->GetSPDS();
  auto fluds =  angleset->fluds;

  const auto num_grps   = angleset->GetNumGrps();
  const auto num_angles = angleset->angles.size();

  //============================== Resize FLUDS non-local incoming Data
  const size_t num_loc_deps = spds.GetLocationDependencies().size();
  if (!upstream_data_initialized)
  {
    angleset->prelocI_outgoing_psi.resize(num_loc_deps, std::vector<double>());
    for (size_t prelocI=0; prelocI<num_loc_deps; prelocI++)
    {
      angleset->prelocI_outgoing_psi[prelocI].resize(
        fluds->prelocI_face_dof_count[prelocI]*num_grps*num_angles,0.0);
    }

    upstream_data_initialized = true;
  }

  //============================== Assume all data is available and now try
  //                               to receive all of it
  bool ready_to_execute = true;
  for (size_t prelocI=0; prelocI<num_loc_deps; prelocI++)
  {
    int locJ = spds.GetLocationDependencies()[prelocI];

    size_t num_mess = prelocI_message_count[prelocI];
    for (int m=0; m<num_mess; m++)
    {
      if (!prelocI_message_received[prelocI][m])
      {
        int message_available = 0;
        MPI_Iprobe(comm_set.MapIonJ(locJ, Chi::mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set.LocICommunicator(Chi::mpi.location_id),
                   &message_available, MPI_STATUS_IGNORE);

        if (not message_available)
        {
          ready_to_execute = false;
          continue;
        }//if message is not available

        //============================ Receive upstream data
        auto& upstream_psi = angleset->prelocI_outgoing_psi[prelocI];

        u_ll_int block_addr   = prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = prelocI_message_size[prelocI][m];

        int error_code =
          MPI_Recv(&upstream_psi[block_addr],
                   static_cast<int>(message_size),
                   MPI_DOUBLE,
                   comm_set.MapIonJ(locJ, Chi::mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set.LocICommunicator(Chi::mpi.location_id),
                   MPI_STATUS_IGNORE);

        prelocI_message_received[prelocI][m] = true;

        if (error_code != MPI_SUCCESS)
        {
          std::stringstream err_stream;
          err_stream << "################# Delayed receive error."
                     << " message size=" << message_size
                     << " as_num=" << angle_set_num
                     << " num_mess=" << num_mess
                     << " m=" << m
                     << " error="
                     << " size=\n";
          char error_string[BUFSIZ];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          MPI_Error_string(error_code, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          Chi::log.LogAllWarning() << err_stream.str();
        }
      }//if not message already received
    }//for message

    if (!ready_to_execute) break;
  }//for predecessor

  if (!ready_to_execute)
    return AngleSetStatus::RECEIVING;
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}