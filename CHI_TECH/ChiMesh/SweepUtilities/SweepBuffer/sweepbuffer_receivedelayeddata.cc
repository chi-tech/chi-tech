#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog     chi_log;
extern ChiMPI     chi_mpi;

//###################################################################
/** Receives delayed data from successor locations. */
void chi_mesh::sweep_management::SweepBuffer::
ReceiveDelayedData(int angle_set_num)
{
  chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();
  MPI_Barrier(MPI_COMM_WORLD);
  //======================================== Receive delayed data
  for (size_t prelocI=0; prelocI<spds->delayed_location_dependencies.size(); prelocI++)
  {
    int locJ = spds->delayed_location_dependencies[prelocI];

    std::vector<double> psi_old(angleset->delayed_prelocI_outgoing_psi[prelocI].size(),0.0);

    for (size_t k=0; k<psi_old.size(); k++)
      psi_old[k] = angleset->delayed_prelocI_outgoing_psi[prelocI][k];

    int num_mess = delayed_prelocI_message_count[prelocI];
    for (int m=0; m<num_mess; m++)
    {

      if (!delayed_prelocI_message_available[prelocI][m])
      {
        int msg_avail = 1;

        MPI_Status status0;
        MPI_Iprobe(comm_set->MapIonJ(locJ,chi_mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi_mpi.location_id],
                   &msg_avail,&status0);

        if (msg_avail != 1)
          chi_log.Log(LOG_ALL)
            << "SweepBuffer: Delayed Data message was not available";

        //============================ Receive upstream data
        u_ll_int block_addr   = delayed_prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = delayed_prelocI_message_size[prelocI][m];

        MPI_Status status;
        int error_code =
          MPI_Recv(&angleset->delayed_prelocI_outgoing_psi[prelocI].data()[block_addr],
                   message_size,
                   MPI_DOUBLE,
                   comm_set->MapIonJ(locJ,chi_mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi_mpi.location_id],
                   &status);

        int num = MPI_Get_count(&status,MPI_DOUBLE,&num);

        if (error_code != MPI_SUCCESS)
        {
          std::stringstream err_stream;
          err_stream << "################# Delayed receive error."
                     << " message size=" << message_size
                     << " as_num=" << angle_set_num
                     << " num_mess=" << num_mess
                     << " m=" << m
                     << " error=" << status.MPI_ERROR
                     << " size=" << num << "\n";
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

    //================================================ Compute norms
    double rel_change = 0.0;
    for (size_t k=0; k<psi_old.size(); k++)
    {
      double psi_new = angleset->delayed_prelocI_outgoing_psi[prelocI][k];
      if (std::fabs(psi_new) >= std::numeric_limits<double>::min())
      {
        double max_change =
          std::fabs(psi_new - psi_old[k])/std::fabs(psi_new);
        rel_change = std::max(max_change,rel_change);
      }
    }
    angleset->delayed_prelocI_norm[prelocI] = rel_change;
  }//for delayed predecessor

  //================================================== Compute norm for local data
  double rel_change = 0.0;
  for (size_t k=0; k<angleset->delayed_local_psi_old.size(); k++)
  {
    double psi_new = angleset->delayed_local_psi[k];
    if (std::fabs(psi_new) >= std::numeric_limits<double>::min())
    {
      double max_change =
        std::fabs(psi_new - angleset->delayed_local_psi_old[k])/std::fabs(psi_new);
      rel_change = std::max(max_change,rel_change);
    }
  }
  angleset->delayed_local_norm = rel_change;
}