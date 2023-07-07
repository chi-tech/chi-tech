#ifndef CHI_ANGLESET_H
#define CHI_ANGLESET_H

#include "mesh/chi_mesh.h"
#include "mesh/SweepUtilities/SweepBuffer/sweepbuffer.h"
#include "mesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_mpi.h>

typedef chi_mesh::sweep_management::SweepBoundary SweepBndry;

#include <memory>

//###################################################################
/**Manages the workstages of a single angle set.*/
class chi_mesh::sweep_management::AngleSet
{
public:
  typedef std::shared_ptr<SweepBndry> SweepBndryPtr;
private:
  size_t                            num_grps;
  const SPDS&                       spds;
  bool                              executed = false;

  chi_mesh::sweep_management::SweepBuffer sweep_buffer;

public:
  std::shared_ptr<FLUDS>             fluds;
  std::vector<size_t>                angles;
  std::map<uint64_t, SweepBndryPtr>& ref_boundaries;
  size_t                             ref_subset;

public:
  AngleSet(size_t in_numgrps,
           size_t in_ref_subset,
           const SPDS& in_spds,
           std::shared_ptr<FLUDS>& in_fluds,
           std::vector<size_t>& angle_indices,
           std::map<uint64_t, std::shared_ptr<SweepBndry>>& sim_boundaries,
           int sweep_eager_limit,
           const chi::ChiMPICommunicatorSet& in_comm_set);

  void InitializeDelayedUpstreamData();

  const chi_mesh::sweep_management::SPDS& GetSPDS() const;

  int GetMaxBufferMessages() const;

  void SetMaxBufferMessages(int new_max);

  size_t GetNumGrps() const;

  AngleSetStatus AngleSetAdvance(
             SweepChunk& sweep_chunk,
             int angle_set_num,
             const std::vector<size_t>& timing_tags,
             ExecutionPermission permission = ExecutionPermission::EXECUTE);
  AngleSetStatus FlushSendBuffers();
  void ResetSweepBuffers();
  bool ReceiveDelayedData(size_t angle_set_num);

  const double* PsiBndry(uint64_t bndry_map,
                         int angle_num,
                         uint64_t cell_local_id,
                         int face_num,
                         int fi,
                         int g,
                         int gs_ss_begin,
                         bool surface_source_active);
  double* ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                     int angle_num,
                                     uint64_t cell_local_id,
                                     int face_num,
                                     int fi,
                                     int gs_ss_begin);
};

#endif //CHI_ANGLESET_H