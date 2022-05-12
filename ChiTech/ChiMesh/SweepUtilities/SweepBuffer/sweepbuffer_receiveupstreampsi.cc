#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog     chi_log;

//###################################################################
/**Check if all upstream dependencies have been met and receives
 * it as it becomes available.*/
chi_mesh::sweep_management::AngleSetStatus
chi_mesh::sweep_management::SweepBuffer::ReceiveUpstreamPsi(int angle_set_num)
{
  auto  spds =  angleset->GetSPDS();
  auto fluds =  angleset->fluds;

  int num_grps   = angleset->GetNumGrps();
  int num_angles = angleset->angles.size();

  //============================== Resize FLUDS non-local incoming Data
  if (!upstream_data_initialized)
  {
    angleset->prelocI_outgoing_psi.resize(
      spds->location_dependencies.size(),std::vector<double>());
    for (size_t prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
    {
      angleset->prelocI_outgoing_psi[prelocI].resize(
        fluds->prelocI_face_dof_count[prelocI]*num_grps*num_angles,0.0);
    }

    upstream_data_initialized = true;
  }

  //============================== Assume all data is available and now try
  //                               to receive all of it
  bool ready_to_execute = true;
  for (size_t prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
  {
    int locJ = spds->location_dependencies[prelocI];

    size_t num_mess = prelocI_message_count[prelocI];
    for (size_t m=0; m<num_mess; m++)
    {

      if (!prelocI_message_available[prelocI][m])
      {
        int msg_avail = 1;

        MPI_Iprobe(comm_set->MapIonJ(locJ,chi::mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi::mpi.location_id],
                   &msg_avail,MPI_STATUS_IGNORE);

        if (msg_avail != 1)
        {
          ready_to_execute = false;
          break;
        }//if message is not available


        //============================ Receive upstream data
        prelocI_message_available[prelocI][m] = true;

        u_ll_int block_addr   = prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = prelocI_message_size[prelocI][m];

        int error_code = MPI_Recv(&angleset->prelocI_outgoing_psi[prelocI].data()[block_addr],
                                  message_size,
                                  MPI_DOUBLE,
                                  comm_set->MapIonJ(locJ,chi::mpi.location_id),
                                  max_num_mess*angle_set_num + m, //tag
                                  comm_set->communicators[chi::mpi.location_id],
                                  MPI_STATUS_IGNORE);

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
          chi_log.Log(LOG_ALLWARNING) << err_stream.str();
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

////###################################################################
///**Check if all upstream dependencies have been met and receives all of
/// it in one go.*/
//chi_mesh::sweep_management::AngleSetStatus
//chi_mesh::sweep_management::SweepBuffer::ReceiveUpstreamPsi(int angle_set_num)
//{
//  auto  spds =  angleset->GetSPDS();
//  auto fluds =  angleset->fluds;
//
//  int num_grps   = angleset->GetNumGrps();
//  int num_angles = angleset->angles.size();
//
//  //============================== Resize FLUDS non-local incoming Data
//  if (!upstream_data_initialized)
//  {
//    angleset->prelocI_outgoing_psi.resize(
//      spds->location_dependencies.size(),std::vector<double>());
//    for (size_t prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
//    {
//      angleset->prelocI_outgoing_psi[prelocI].resize(
//        fluds->prelocI_face_dof_count[prelocI]*num_grps*num_angles,0.0);
//    }
//
//    upstream_data_initialized = true;
//  }
//
//  //============================== Assume all data is available and now try
//  //                               to disprove it
//  bool ready_to_execute = true;
//  for (size_t prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
//  {
//    int locJ = spds->location_dependencies[prelocI];
//
//    size_t num_mess = prelocI_message_count[prelocI];
//    for (size_t m=0; m<num_mess; m++)
//    {
//
//      if (!prelocI_message_available[prelocI][m])
//      {
//        int msg_avail = 1;
//
//        MPI_Iprobe(comm_set->MapIonJ(locJ,chi::mpi.location_id),
//                   max_num_mess*angle_set_num + m, //tag
//                   comm_set->communicators[chi::mpi.location_id],
//                   &msg_avail,MPI_STATUS_IGNORE);
//
//        if (msg_avail != 1)
//        {
//          ready_to_execute = false;
//          break;
//        }//if message is not available
//      }//if not message already received
//    }//for message
//
//    if (!ready_to_execute) break;
//  }//for predecessor
//
//  //============================== If all data is available, receive it
//  if (ready_to_execute)
//  {
//    for (size_t prelocI=0; prelocI<spds->location_dependencies.size(); prelocI++)
//    {
//      int locJ = spds->location_dependencies[prelocI];
//
//      size_t num_mess = prelocI_message_count[prelocI];
//      for (size_t m=0; m<num_mess; m++)
//      {
//        //Prevent rereceive
//        if (prelocI_message_available[prelocI][m]) continue;
//        prelocI_message_available[prelocI][m] = true;
//
//        u_ll_int block_addr   = prelocI_message_blockpos[prelocI][m];
//        u_ll_int message_size = prelocI_message_size[prelocI][m];
//
//        int error_code = MPI_Recv(&angleset->prelocI_outgoing_psi[prelocI].data()[block_addr],
//                                  message_size,
//                                  MPI_DOUBLE,
//                                  comm_set->MapIonJ(locJ,chi::mpi.location_id),
//                                  max_num_mess*angle_set_num + m, //tag
//                                  comm_set->communicators[chi::mpi.location_id],
//                                  MPI_STATUS_IGNORE);
//
//        if (error_code != MPI_SUCCESS)
//        {
//          std::stringstream err_stream;
//          err_stream << "################# Delayed receive error."
//                     << " message size=" << message_size
//                     << " as_num=" << angle_set_num
//                     << " num_mess=" << num_mess
//                     << " m=" << m
//                     << " error="
//                     << " size=\n";
//          char error_string[BUFSIZ];
//          int length_of_error_string, error_class;
//          MPI_Error_class(error_code, &error_class);
//          MPI_Error_string(error_class, error_string, &length_of_error_string);
//          err_stream << error_string << "\n";
//          MPI_Error_string(error_code, error_string, &length_of_error_string);
//          err_stream << error_string << "\n";
//          chi_log.Log(LOG_ALLWARNING) << err_stream.str();
//        }
//      }//for message
//
//    }//for predecessor
//  }
//
//
//  if (!ready_to_execute)
//    return AngleSetStatus::RECEIVING;
//  else
//    return AngleSetStatus::READY_TO_EXECUTE;
//}


