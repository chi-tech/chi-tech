#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/** Receives delayed data from successor locations. */
void chi_mesh::sweep_management::SweepBuffer::
  ReceiveDelayedData(int angle_set_num)
{
  auto spds =  angleset->GetSPDS();

  //======================================== Receive delayed data
  const size_t num_delayed_loc_deps = spds->delayed_location_dependencies.size();
  for (size_t prelocI=0; prelocI<num_delayed_loc_deps; prelocI++)
  {
    int locJ = spds->delayed_location_dependencies[prelocI];

    auto& outgoing_psi = angleset->delayed_prelocI_outgoing_psi[prelocI];

    const size_t psi_size = outgoing_psi.size();
    std::vector<double> psi_old(psi_size,0.0);

    for (size_t k=0; k<psi_size; k++)
      psi_old[k] = outgoing_psi[k];

    int num_mess = delayed_prelocI_message_count[prelocI];
    for (int m=0; m<num_mess; m++)
    {
      MPI_Probe(comm_set->MapIonJ(locJ,chi::mpi.location_id),
                max_num_mess*angle_set_num + m, //tag
                comm_set->communicators[chi::mpi.location_id],
                nullptr);

      //============================ Receive upstream data
      u_ll_int block_addr   = delayed_prelocI_message_blockpos[prelocI][m];
      u_ll_int message_size = delayed_prelocI_message_size[prelocI][m];

      MPI_Status status;
      int error_code =
        MPI_Recv(&outgoing_psi[block_addr],
                 static_cast<int>(message_size),
                 MPI_DOUBLE,
                 comm_set->MapIonJ(locJ,chi::mpi.location_id),
                 max_num_mess*angle_set_num + m, //tag
                 comm_set->communicators[chi::mpi.location_id],
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
        chi::log.LogAllWarning() << err_stream.str();
      }
    }//for message
  }//for delayed predecessor

  //================================================== Copy non-local delayed Psi
  //                                                   to Psi_old
  for (size_t prelocI=0; prelocI<num_delayed_loc_deps; prelocI++)
  {
    auto& psi_new = angleset->delayed_prelocI_outgoing_psi[prelocI];
    auto& psi_old = angleset->delayed_prelocI_outgoing_psi_old[prelocI];
    psi_old = psi_new;
  }

  //================================================== Copy local delayed Psi
  //                                                   to Psi_old
  angleset->delayed_local_psi_old = angleset->delayed_local_psi;

}