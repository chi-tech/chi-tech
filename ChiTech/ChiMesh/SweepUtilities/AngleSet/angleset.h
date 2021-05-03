#ifndef CHI_ANGLESET_H
#define CHI_ANGLESET_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/SweepUtilities/SweepBuffer/sweepbuffer.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_mpi.h>

typedef chi_mesh::sweep_management::BoundaryBase SweepBndry;

#include <memory>

//###################################################################
/**Manages the workstages of a single angle set.*/
class chi_mesh::sweep_management::AngleSet
{
private:
  int                               num_grps;
  std::shared_ptr<SPDS>             spds;
  bool                              executed;

  chi_mesh::sweep_management::SweepBuffer sweep_buffer;

public:
  FLUDS*                            fluds;
  std::vector<int>                  angles;
  std::vector<std::shared_ptr<SweepBndry>>&         ref_boundaries;
  int                               ref_subset;

  //FLUDS
  std::vector<std::vector<double>>  local_psi;
  std::vector<double>               delayed_local_psi;
  std::vector<double>               delayed_local_psi_old;
  std::vector<std::vector<double>>  deplocI_outgoing_psi;
  std::vector<std::vector<double>>  prelocI_outgoing_psi;
  std::vector<std::vector<double>>  boundryI_incoming_psi;

  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi;
  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi_old;
  std::vector<double>               delayed_prelocI_norm;
  double                            delayed_local_norm;

  AngleSet(int in_numgrps,
           int in_ref_subset,
           std::shared_ptr<SPDS>& in_spds,
           FLUDS* in_fluds,
           std::vector<int>& angle_indices,
           std::vector<std::shared_ptr<SweepBndry>>& sim_boundaries,
           int sweep_eager_limit,
           ChiMPICommunicatorSet* in_comm_set);

  void InitializeDelayedUpstreamData();

  std::shared_ptr<chi_mesh::sweep_management::SPDS> GetSPDS();

  int GetMaxBufferMessages();

  void SetMaxBufferMessages(int new_max);

  int GetNumGrps();

  AngleSetStatus AngleSetAdvance(
             SweepChunk *sweep_chunk,
             int angle_set_num,
             const std::vector<size_t>& timing_tags,
             ExecutionPermission permission = ExecutionPermission::EXECUTE);
  AngleSetStatus FlushSendBuffers();
  void ResetSweepBuffers();
  void ReceiveDelayedData(int angle_set_num);

  double* PsiBndry(int bndry_map,
                   int angle_num,
                   int cell_local_id,
                   int face_num,
                   int fi,
                   int g,
                   int gs_ss_begin,
                   bool suppress_surface_src);
  double* ReflectingPsiOutBoundBndry(int bndry_map,
                                     int angle_num,
                                     int cell_local_id,
                                     int face_num,
                                     int fi,
                                     int gs_ss_begin);
};

#endif //CHI_ANGLESET_H