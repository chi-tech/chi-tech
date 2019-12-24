#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <chi_mpi.h>
extern ChiMPI chi_mpi;

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**AngleSet constructor.*/
chi_mesh::sweep_management::AngleSet::
AngleSet(int in_numgrps,
         int in_ref_subset,
         SPDS* in_spds,
         std::vector<int>& angle_indices,
         std::vector<SweepBndry*>& sim_boundaries,
         int sweep_eager_limit,
         ChiMPICommunicatorSet* in_comm_set):
  sweep_buffer(this,sweep_eager_limit,in_comm_set),
  ref_boundaries(sim_boundaries)
{
  num_grps = in_numgrps;
  spds     = in_spds;
  executed = false;
  ref_subset = in_ref_subset;
  std::copy(angle_indices.begin(),
            angle_indices.end(),
            std::back_inserter(angles));

  auto primary_fluds = new chi_mesh::sweep_management::PRIMARY_FLUDS(num_grps);
  primary_fluds->InitializeAlphaElements(spds);
  primary_fluds->InitializeBetaElements(spds);

  fluds = primary_fluds;

  sweep_buffer.BuildMessageStructure();

  delayed_local_norm = 0.0;
};

//###################################################################
/**AngleSet constructor.*/
chi_mesh::sweep_management::AngleSet::
AngleSet(int in_numgrps,
         int in_ref_subset,
         SPDS* in_spds,
         FLUDS* in_fluds,
         std::vector<int>& angle_indices,
         std::vector<SweepBndry*>& sim_boundaries,
         int sweep_eager_limit,
         ChiMPICommunicatorSet* in_comm_set):
  sweep_buffer(this,sweep_eager_limit,in_comm_set),
  ref_boundaries(sim_boundaries)
{
  num_grps = in_numgrps;
  spds     = in_spds;
  executed = false;
  ref_subset = in_ref_subset;
  std::copy(angle_indices.begin(),
            angle_indices.end(),
            std::back_inserter(angles));

  fluds = in_fluds;

  sweep_buffer.BuildMessageStructure();

  delayed_local_norm = 0.0;
};

//###################################################################
/**Initializes delayed upstream data. This method gets called
 * when a sweep scheduler is constructed.*/
void chi_mesh::sweep_management::AngleSet::
InitializeDelayedUpstreamData()
{
  sweep_buffer.InitializeDelayedUpstreamData();
}

//###################################################################
/**This function advances the work stages of an angleset.*/
chi_mesh::sweep_management::AngleSetStatus
chi_mesh::sweep_management::AngleSet::
AngleSetAdvance(chi_mesh::sweep_management::SweepChunk *sweep_chunk,
                int angle_set_num,
                const std::vector<size_t>& timing_tags,
                chi_mesh::sweep_management::ExecutionPermission permission)
{
  typedef AngleSetStatus Status;

  if (executed)
  {
    if (!sweep_buffer.DoneSending())
      sweep_buffer.ClearDownstreamBuffers();
    return AngleSetStatus::FINISHED;
  }

  //Check upstream data available
  Status status = sweep_buffer.ReceiveUpstreamPsi(angle_set_num);

  //Also check boundaries
  for (auto bndry : ref_boundaries)
    if (not bndry->CheckAnglesReadyStatus(angles,ref_subset))
      {status = Status::RECEIVING; break;}

  if      (status == Status::RECEIVING) return status;
  else if (status == Status::READY_TO_EXECUTE and
           permission == ExecutionPermission::EXECUTE)
  {
    sweep_buffer.InitializeLocalAndDownstreamBuffers();

    chi_log.LogEvent(timing_tags[0],ChiLog::EventType::EVENT_BEGIN);
    sweep_chunk->Sweep(this); //Execute chunk
    chi_log.LogEvent(timing_tags[0],ChiLog::EventType::EVENT_END);

    //Send outgoing psi and clear local and receive buffers
    sweep_buffer.SendDownstreamPsi(angle_set_num);
    sweep_buffer.ClearLocalAndReceiveBuffers();

    //Update boundary readiness
    for (auto bndry : ref_boundaries)
      bndry->UpdateAnglesReadyStatus(angles,ref_subset);

    executed = true;
    return AngleSetStatus::FINISHED;
  }
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}

//###################################################################
/**Returns a reference to the associated spds.*/
chi_mesh::sweep_management::SPDS*
  chi_mesh::sweep_management::AngleSet::GetSPDS()
{
  return spds;
}

//###################################################################
/**Returns the maximum buffer size from the sweepbuffer.*/
int chi_mesh::sweep_management::AngleSet::GetMaxBufferMessages()
{
  return sweep_buffer.max_num_mess;
}

//###################################################################
/**Sets the maximum buffer size for the sweepbuffer.*/
void chi_mesh::sweep_management::AngleSet::SetMaxBufferMessages(int new_max)
{
  sweep_buffer.max_num_mess = new_max;
}

//###################################################################
/**Returns the number of groups associated with the angleset.*/
int chi_mesh::sweep_management::AngleSet::GetNumGrps()
{
  return num_grps;
}

//###################################################################
/**Resets the sweep buffer.*/
void chi_mesh::sweep_management::AngleSet::ResetSweepBuffers()
{
  sweep_buffer.Reset();
  executed = false;
}

//###################################################################
/**Instructs the sweep buffer to receive delayed data.*/
void chi_mesh::sweep_management::AngleSet::ReceiveDelayedData(int angle_set_num)
{
  sweep_buffer.ReceiveDelayedData(angle_set_num);
}

//###################################################################
/**Returns a pointer to a boundary flux data.*/
double* chi_mesh::sweep_management::AngleSet::
PsiBndry(int bndry_map,
         int angle_num,
         int cell_local_id,
         int face_num,
         int fi,
         int g,
         int gs_ss_begin,
         bool suppress_surface_src)
{
  double* Psi = &ref_boundaries[bndry_map]->boundary_flux[g];

  if (suppress_surface_src)
    Psi = &ref_boundaries[bndry_map]->zero_boundary_flux[g];

  if (ref_boundaries[bndry_map]->IsReflecting())
  {
    Psi = ref_boundaries[bndry_map]->HeterogenousPsiIncoming(
      angle_num, cell_local_id, face_num, fi, gs_ss_begin);

  }

  return Psi;
}

//###################################################################
/**Returns a pointer to outbound boundary flux data.*/
double* chi_mesh::sweep_management::AngleSet::
ReflectingPsiOutBoundBndry(int bndry_map,
                           int angle_num,
                           int cell_local_id,
                           int face_num,
                           int fi,
                           int gs_ss_begin)
{
  return ref_boundaries[bndry_map]->HeterogenousPsiOutgoing(angle_num,
                                                            cell_local_id,
                                                            face_num,
                                                            fi,
                                                            gs_ss_begin);
}