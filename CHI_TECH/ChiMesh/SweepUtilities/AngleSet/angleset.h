#ifndef _chi_angleset_h
#define _chi_angleset_h

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/SweepUtilities/SweepBuffer/sweepbuffer.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include <ChiMesh/SweepUtilities/FLUDS/FLUDS.h>

#include <chi_mpi.h>

typedef chi_mesh::sweep_management::BoundaryBase SweepBndry;

//###################################################################
/**Manages the workstages of a single angle set.*/
class chi_mesh::sweep_management::AngleSet
{
private:
  int              num_grps;
  SPDS*            spds;
  bool             executed;

  chi_mesh::sweep_management::SweepBuffer sweep_buffer;

public:
  FLUDS*                                 fluds;
  std::vector<int>                       angles;
  std::vector<SweepBndry*>&              ref_boundaries;
  int                                    ref_subset;

  //FLUDS
  std::vector<std::vector<double>>  local_psi;
  std::vector<double>               delayed_local_psi;
  std::vector<double>               delayed_local_psi_old;
  std::vector<std::vector<double>>  deplocI_outgoing_psi;
  std::vector<std::vector<double>>  prelocI_outgoing_psi;
  std::vector<std::vector<double>>  boundryI_incoming_psi;

  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi;
  std::vector<double>               delayed_prelocI_norm;
  double                            delayed_local_norm;

  AngleSet(int in_numgrps,
           int in_ref_subset,
           SPDS* in_spds,
           std::vector<int>& angle_indices,
           std::vector<SweepBndry*>& sim_boundaries,
           int sweep_eager_limit,
           ChiMPICommunicatorSet* in_comm_set);

  AngleSet(int in_numgrps,
           int in_ref_subset,
           SPDS* in_spds,
           FLUDS* in_fluds,
           std::vector<int>& angle_indices,
           std::vector<SweepBndry*>& sim_boundaries,
           int sweep_eager_limit,
           ChiMPICommunicatorSet* in_comm_set);

  void InitializeDelayedUpstreamData();

  SPDS* GetSPDS();

  int GetMaxBufferMessages();

  void SetMaxBufferMessages(int new_max);

  int GetNumGrps();

  AngleSetStatus AngleSetAdvance(
             SweepChunk *sweep_chunk,
             int angle_set_num,
             const std::vector<size_t>& timing_tags,
             ExecutionPermission permission = ExecutionPermission::EXECUTE);
  void ResetSweepBuffers();
  void ReceiveDelayedData(int angle_set_num);

  double* PsiBndry(int bndry_face_count, int bndry_map,
                      int face_dof, int g,int angle_num);

};

#endif