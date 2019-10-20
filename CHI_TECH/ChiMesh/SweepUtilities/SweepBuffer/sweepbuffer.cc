#include "sweepbuffer.h"

#include "ChiMesh/SweepUtilities/AngleSet/angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_log.h>
#include <chi_mpi.h>
#include <ChiConsole/chi_console.h>

extern ChiLog     chi_log;
extern ChiMPI     chi_mpi;
extern ChiConsole chi_console;

extern double chi_global_timings[20];

//###################################################################
/**Constructor.*/
chi_mesh::sweep_management::SweepBuffer::
  SweepBuffer(chi_mesh::sweep_management::AngleSet* ref_angleset,
              int sweep_eager_limit,
              ChiMPICommunicatorSet* in_comm_set)
{
  angleset = ref_angleset;
  initialized = false;
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;
  EAGER_LIMIT = sweep_eager_limit;
  comm_set = in_comm_set;

  max_num_mess = 0;
}

//###################################################################
/**Builds message structure.
 *
 * Outgoing and incoming data needs to be sub-divided into messages
 * each of which is smaller than the MPI eager-limit. There are
 * three parts to this: predecessors, delayed-predecessors and successors.*/
void chi_mesh::sweep_management::SweepBuffer::BuildMessageStructure()
{
//============================================= Check angleset is complete
  if (angleset->angles.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "A call to SweepBuffer::BuildMessageStructure() has been made without"
         " an initialized angleset.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();
  chi_mesh::sweep_management::FLUDS* fluds=  angleset->fluds;

  int num_grps   = angleset->GetNumGrps();
  int num_angles = angleset->angles.size();

  //============================================= Predecessor locations
  size_t num_dependencies = spds->location_dependencies.size();

  prelocI_message_count.resize(num_dependencies,0);
  prelocI_message_size.resize(num_dependencies);
  prelocI_message_blockpos.resize(num_dependencies);
  prelocI_message_available.clear();

  for (size_t prelocI=0; prelocI<num_dependencies; prelocI++)
  {
    u_ll_int num_unknowns =
      fluds->prelocI_face_dof_count[prelocI]*num_grps*num_angles;

    u_ll_int message_size  = num_unknowns;
    int      message_count = 1;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = num_angles;
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns*8/(double)(double)EAGER_LIMIT);
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }

    prelocI_message_count[prelocI] = message_count;

    u_ll_int pre_block_pos = 0;
    for (int m=0; m<(message_count-1); m++)
    {
      prelocI_message_size[prelocI].push_back(message_size);
      prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns>0)
    {
      prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      prelocI_message_size[prelocI].push_back(num_unknowns);
    }

    prelocI_message_available.emplace_back(message_count,false);
  }

  //============================================= Delayed Predecessor locations
  size_t num_delayed_dependencies = spds->delayed_location_dependencies.size();

  delayed_prelocI_message_count.resize(num_delayed_dependencies,0);
  delayed_prelocI_message_size.resize(num_delayed_dependencies);
  delayed_prelocI_message_blockpos.resize(num_delayed_dependencies);
  delayed_prelocI_message_available.clear();

  for (size_t prelocI=0; prelocI<num_delayed_dependencies; prelocI++)
  {
    angleset->delayed_prelocI_norm.push_back(0.0);

    u_ll_int num_unknowns =
      fluds->delayed_prelocI_face_dof_count[prelocI]*num_grps*num_angles;

    u_ll_int message_size  = num_unknowns;
    int      message_count = 1;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = num_angles;
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns*8/(double)(double)EAGER_LIMIT);
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }

    delayed_prelocI_message_count[prelocI] = message_count;

    u_ll_int pre_block_pos = 0;
    for (int m=0; m<(message_count-1); m++)
    {
      delayed_prelocI_message_size[prelocI].push_back(message_size);
      delayed_prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns>0)
    {
      delayed_prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      delayed_prelocI_message_size[prelocI].push_back(num_unknowns);
    }

    delayed_prelocI_message_available.emplace_back(message_count,false);
  }


  //============================================= Successor locations
  size_t num_successors = spds->location_successors.size();

  deplocI_message_count.resize(num_successors,0);
  deplocI_message_size.resize(num_successors);
  deplocI_message_blockpos.resize(num_successors);

  deplocI_message_sent.clear();
  deplocI_message_request.clear();

  for (size_t deplocI=0; deplocI<num_successors; deplocI++)
  {
    u_ll_int num_unknowns =
      fluds->deplocI_face_dof_count[deplocI]*num_grps*num_angles;

    u_ll_int message_size  = num_unknowns;
    int      message_count = 1;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = num_angles;
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns*8/(double)(double)EAGER_LIMIT);
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }

    chi_global_timings[12] += num_unknowns;
    chi_global_timings[13] += 1.0;
    chi_global_timings[14] += message_count;
    chi_global_timings[15] += 1.0;

    deplocI_message_count[deplocI] = message_count;

    u_ll_int dep_block_pos = 0;
    for (int m=0; m<(message_count-1); m++)
    {
      deplocI_message_size[deplocI].push_back(message_size);
      deplocI_message_blockpos[deplocI].push_back(dep_block_pos);
      num_unknowns -= message_size;
      dep_block_pos += message_size;
    }
    if (num_unknowns>0)
    {
      deplocI_message_blockpos[deplocI].push_back(dep_block_pos);
      deplocI_message_size[deplocI].push_back(num_unknowns);
    }

    deplocI_message_sent.emplace_back(message_count,false);
    deplocI_message_request.emplace_back(message_count,MPI_Request());
  }

  angleset->fluds->SetReferencePsi(&angleset->local_psi,
                                   &angleset->delayed_local_psi,
                                   &angleset->deplocI_outgoing_psi,
                                   &angleset->prelocI_outgoing_psi,
                                   &angleset->boundryI_incoming_psi,
                                   &angleset->delayed_prelocI_outgoing_psi);

  //================================================== All reduce to get
  //                                                   maximum message count
  int angset_max_message_count = 0;
  for (size_t prelocI=0; prelocI<num_dependencies; prelocI++)
    angset_max_message_count =
      std::max(prelocI_message_count[prelocI], angset_max_message_count);

  for (size_t prelocI=0; prelocI<num_delayed_dependencies; prelocI++)
    angset_max_message_count =
      std::max(delayed_prelocI_message_count[prelocI], angset_max_message_count);

  for (size_t deplocI=0; deplocI<num_successors; deplocI++)
    angset_max_message_count =
      std::max(deplocI_message_count[deplocI], angset_max_message_count);

  //Temporarily assign max_num_mess tot he local maximum
  max_num_mess = angset_max_message_count;
}


//###################################################################
/** This is the final level of initialization before a sweep-chunk executes.
 * Once all upstream dependencies are met and if the sweep scheduler places
 * this angleset as "ready-to-execute", then the angle-set will call this
 * method. It is also fairly important in terms of memory to only allocate
 * these chunks of memory since they form the majority of memory requirements.*/
void chi_mesh::sweep_management::SweepBuffer::
  InitializeBuffers()
{
  if (!data_initialized)
  {
    chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();
    chi_mesh::sweep_management::FLUDS* fluds=  angleset->fluds;

    int num_grps   = angleset->GetNumGrps();
    int num_angles = angleset->angles.size();

    //============================ Resize FLUDS local outgoing Data
    angleset->local_psi.resize(fluds->num_face_categories);
    // fc = face category
    for (size_t fc = 0; fc<fluds->num_face_categories; fc++)
    {
      angleset->local_psi[fc].resize(fluds->local_psi_stride[fc]*
                                     fluds->local_psi_max_elements[fc]*
                                     num_grps*num_angles,0.0);
    }

    //============================ Resize FLUDS non-local outgoing Data
    angleset->deplocI_outgoing_psi.resize(
      spds->location_successors.size(),std::vector<double>());
    for (size_t deplocI=0; deplocI<spds->location_successors.size(); deplocI++)
    {
      angleset->deplocI_outgoing_psi[deplocI].resize(
        fluds->deplocI_face_dof_count[deplocI]*num_grps*num_angles,0.0);
    }

    //================================================ Make a memory query
    double memory_mb = chi_console.GetMemoryUsageInMB();

    if (memory_mb>chi_global_timings[9])
      chi_global_timings[9] = memory_mb;

    chi_global_timings[10] += memory_mb;
    chi_global_timings[11] += 1.0;

    data_initialized = true;
  }

  //================================================== Copy delayed Psi to Psi_old
  angleset->delayed_local_psi_old.clear();
  std::copy(angleset->delayed_local_psi.begin(),
            angleset->delayed_local_psi.end(),
            std::back_inserter(angleset->delayed_local_psi_old));
}

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

//###################################################################
/**Sends downstream psi.*/
void chi_mesh::sweep_management::SweepBuffer::
 ClearDownstreamBuffers()
{
  if (done_sending) return;

  chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();

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
/**Check if all upstream dependencies have been met.*/
chi_mesh::sweep_management::AngleSetStatus
chi_mesh::sweep_management::SweepBuffer::ReceiveUpstreamPsi(int angle_set_num)
{
  chi_mesh::sweep_management::SPDS*  spds =  angleset->GetSPDS();
  chi_mesh::sweep_management::FLUDS* fluds=  angleset->fluds;

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

        MPI_Iprobe(comm_set->MapIonJ(locJ,chi_mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi_mpi.location_id],
                   &msg_avail,MPI_STATUS_IGNORE);

        if (msg_avail != 1)
        {
          ready_to_execute = false;
          break;
        }//if message is not available
        //============================ Receive upstream data
        else
        {
          prelocI_message_available[prelocI][m] = true;

          u_ll_int block_addr   = prelocI_message_blockpos[prelocI][m];
          u_ll_int message_size = prelocI_message_size[prelocI][m];

          int error_code = MPI_Recv(&angleset->prelocI_outgoing_psi[prelocI].data()[block_addr],
                   message_size,
                   MPI_DOUBLE,
                   comm_set->MapIonJ(locJ,chi_mpi.location_id),
                   max_num_mess*angle_set_num + m, //tag
                   comm_set->communicators[chi_mpi.location_id],
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


//###################################################################
/**Receive all upstream Psi.*/
void chi_mesh::sweep_management::SweepBuffer::
  ClearLocalAndReceiveBuffers()
{
  auto empty_vector = std::vector<std::vector<double>>(0);
  angleset->local_psi.swap(empty_vector);

  empty_vector = std::vector<std::vector<double>>(0);
  angleset->prelocI_outgoing_psi.swap(empty_vector);

//  angleset->local_psi.clear();
//  angleset->local_psi.shrink_to_fit();
//
//  //Clear receive buffers
//  angleset->prelocI_outgoing_psi.clear();
//  angleset->prelocI_outgoing_psi.shrink_to_fit();
}