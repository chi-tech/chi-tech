#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "angleset.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/sweepchunk_base.h"

#include "chi_mpi.h"
#include "chi_log.h"


//###################################################################
/**AngleSet constructor.*/
chi_mesh::sweep_management::AngleSet::
AngleSet(size_t in_numgrps,
         size_t in_ref_subset,
         const SPDS& in_spds,
         FLUDS* in_fluds,
         std::vector<size_t>& angle_indices,
         std::map<uint64_t, std::shared_ptr<SweepBndry>>& sim_boundaries,
         int sweep_eager_limit,
         const chi_objects::ChiMPICommunicatorSet& in_comm_set):
  num_grps(in_numgrps),
  spds(in_spds),
  sweep_buffer(this,sweep_eager_limit,in_comm_set),
  fluds(in_fluds),
  angles(angle_indices),
  ref_boundaries(sim_boundaries),
  ref_subset(in_ref_subset)
{
  sweep_buffer.BuildMessageStructure();
}

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
  AngleSetAdvance(chi_mesh::sweep_management::SweepChunk& sweep_chunk,
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
  for (auto& [bid,bndry] : ref_boundaries)
    if (not bndry->CheckAnglesReadyStatus(angles,ref_subset))
      {status = Status::RECEIVING; break;}

  if      (status == Status::RECEIVING) return status;
  else if (status == Status::READY_TO_EXECUTE and
           permission == ExecutionPermission::EXECUTE)
  {
    sweep_buffer.InitializeLocalAndDownstreamBuffers();

    chi::log.LogEvent(timing_tags[0], chi_objects::ChiLog::EventType::EVENT_BEGIN);
    sweep_chunk.Sweep(this); //Execute chunk
    chi::log.LogEvent(timing_tags[0], chi_objects::ChiLog::EventType::EVENT_END);

    //Send outgoing psi and clear local and receive buffers
    sweep_buffer.SendDownstreamPsi(angle_set_num);
    sweep_buffer.ClearLocalAndReceiveBuffers();

    //Update boundary readiness
    for (auto& [bid,bndry] : ref_boundaries)
      bndry->UpdateAnglesReadyStatus(angles,ref_subset);

    executed = true;
    return AngleSetStatus::FINISHED;
  }
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}

//###################################################################
/***/
chi_mesh::sweep_management::AngleSetStatus
  chi_mesh::sweep_management::AngleSet::FlushSendBuffers()
{
  if (!sweep_buffer.DoneSending())
    sweep_buffer.ClearDownstreamBuffers();

  if (sweep_buffer.DoneSending())
    return AngleSetStatus::MESSAGES_SENT;

  return AngleSetStatus::MESSAGES_PENDING;
}

//###################################################################
/**Returns a reference to the associated spds.*/
const chi_mesh::sweep_management::SPDS&
  chi_mesh::sweep_management::AngleSet::GetSPDS() const
{
  return spds;
}

//###################################################################
/**Returns the maximum buffer size from the sweepbuffer.*/
int chi_mesh::sweep_management::AngleSet::GetMaxBufferMessages() const
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
size_t chi_mesh::sweep_management::AngleSet::GetNumGrps() const
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
bool chi_mesh::sweep_management::AngleSet::ReceiveDelayedData(size_t angle_set_num)
{
  return sweep_buffer.ReceiveDelayedData(static_cast<int>(angle_set_num));
}

//###################################################################
/**Returns a pointer to a boundary flux data.*/
const double* chi_mesh::sweep_management::AngleSet::
PsiBndry(uint64_t bndry_map,
         int angle_num,
         uint64_t cell_local_id,
         int face_num,
         int fi,
         int g,
         int gs_ss_begin,
         bool surface_source_active)
{
  if (ref_boundaries[bndry_map]->IsReflecting())
    return ref_boundaries[bndry_map]->HeterogeneousPsiIncoming(
        cell_local_id, face_num, fi, angle_num, g, gs_ss_begin);

  if (not surface_source_active)
    return ref_boundaries[bndry_map]->ZeroFlux(g);

  return ref_boundaries[bndry_map]->HeterogeneousPsiIncoming(
      cell_local_id, face_num, fi, angle_num, g, gs_ss_begin);
}

//###################################################################
/**Returns a pointer to outbound boundary flux data.*/
double* chi_mesh::sweep_management::AngleSet::
ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                           int angle_num,
                           uint64_t cell_local_id,
                           int face_num,
                           int fi,
                           int gs_ss_begin)
{
  return ref_boundaries[bndry_map]->HeterogeneousPsiOutgoing(cell_local_id,
                                                             face_num,
                                                             fi,
                                                             angle_num,
                                                             gs_ss_begin);
}