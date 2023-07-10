#include "sweepbuffer.h"

#include "mesh/SweepUtilities/AngleSet/angleset.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
#include "console/chi_console.h"

//###################################################################
/**Builds message structure.
 *
 * Outgoing and incoming data needs to be sub-divided into messages
 * each of which is smaller than the MPI eager-limit. There are
 * three parts to this: predecessors, delayed-predecessors and successors.
 *
 * This method gets called by an angleset that subscribes to this
 * sweepbuffer.*/
void chi_mesh::sweep_management::SweepBuffer::BuildMessageStructure()
{
//============================================= Check angleset is complete
  if (angleset->angles.empty())
  {
    Chi::log.LogAllError()
      << "A call to SweepBuffer::BuildMessageStructure() has been made without"
         " an initialized angleset.";
    Chi::Exit(EXIT_FAILURE);
  }

  const auto& spds =  angleset->GetSPDS();
  auto fluds=  angleset->fluds;

  const auto num_grps   = angleset->GetNumGrps();
  const auto num_angles = angleset->angles.size();

  //============================================= Predecessor locations
  size_t num_dependencies = spds.location_dependencies.size();

  prelocI_message_count.resize(num_dependencies,0);
  prelocI_message_size.resize(num_dependencies);
  prelocI_message_blockpos.resize(num_dependencies);
  prelocI_message_received.clear();

  for (size_t prelocI=0; prelocI<num_dependencies; prelocI++)
  {
    u_ll_int num_unknowns =
      fluds->prelocI_face_dof_count[prelocI]*num_grps*num_angles;

    u_ll_int message_size;
    int      message_count;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles);
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

    prelocI_message_received.emplace_back(message_count, false);
  }//for prelocI

  //============================================= Delayed Predecessor locations
  size_t num_delayed_dependencies = spds.delayed_location_dependencies.size();

  delayed_prelocI_message_count.resize(num_delayed_dependencies,0);
  delayed_prelocI_message_size.resize(num_delayed_dependencies);
  delayed_prelocI_message_blockpos.resize(num_delayed_dependencies);
  delayed_prelocI_message_received.clear();

  for (size_t prelocI=0; prelocI<num_delayed_dependencies; prelocI++)
  {
    u_ll_int num_unknowns =
      fluds->delayed_prelocI_face_dof_count[prelocI]*num_grps*num_angles;

    u_ll_int message_size;
    int      message_count;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles);
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

    delayed_prelocI_message_received.emplace_back(message_count, false);
  }


  //============================================= Successor locations
  size_t num_successors = spds.location_successors.size();

  deplocI_message_count.resize(num_successors,0);
  deplocI_message_size.resize(num_successors);
  deplocI_message_blockpos.resize(num_successors);

  deplocI_message_request.clear();

  for (size_t deplocI=0; deplocI<num_successors; deplocI++)
  {
    u_ll_int num_unknowns =
      fluds->deplocI_face_dof_count[deplocI]*num_grps*num_angles;

    u_ll_int message_size;
    int      message_count;
    if ((num_unknowns*8)<=EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles);
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns*8/(double)(double)EAGER_LIMIT);
      message_size  = ceil((double)num_unknowns/(double)message_count);
    }

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

    deplocI_message_request.emplace_back(message_count,MPI_Request());
  }

  angleset->fluds->SetReferencePsi(&angleset->local_psi,
                                   &angleset->delayed_local_psi,
                                   &angleset->delayed_local_psi_old,
                                   &angleset->deplocI_outgoing_psi,
                                   &angleset->prelocI_outgoing_psi,
                                   &angleset->boundryI_incoming_psi,
                                   &angleset->delayed_prelocI_outgoing_psi,
                                   &angleset->delayed_prelocI_outgoing_psi_old);

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