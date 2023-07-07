#ifndef CHI_ANGLESET_H
#define CHI_ANGLESET_H

#include "mesh/chi_mesh.h"
#include "mesh/SweepUtilities/Communicators/AAH_AsynComm.h"
#include "mesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include "mesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_mpi.h>

typedef chi_mesh::sweep_management::SweepBoundary SweepBndry;

#include <memory>

namespace chi_mesh::sweep_management
{

class AngleSet
{
public:
  typedef std::shared_ptr<SweepBndry> SweepBndryPtr;

  AngleSet(size_t num_groups, const SPDS& spds,
           std::shared_ptr<FLUDS>& fluds,
           const std::vector<size_t>& angle_indices,
           std::map<uint64_t, SweepBndryPtr>& sim_boundaries,
           size_t in_ref_subset);


  const SPDS& GetSPDS() const;
  FLUDS& GetFLUDS();
  size_t GetRefGroupSubset() const;
  const std::vector<size_t>& GetAngleIndices() const;
  std::map<uint64_t, SweepBndryPtr>& GetBoundaries();

  // Virtual methods
  virtual void InitializeDelayedUpstreamData() = 0;

  virtual int GetMaxBufferMessages() const = 0;

  virtual void SetMaxBufferMessages(int new_max) = 0;

  virtual AngleSetStatus AngleSetAdvance(
    SweepChunk& sweep_chunk,
    int angle_set_num,
    const std::vector<size_t>& timing_tags,
    ExecutionPermission permission) = 0;
  virtual AngleSetStatus FlushSendBuffers() = 0;
  virtual void ResetSweepBuffers() = 0;
  virtual bool ReceiveDelayedData(size_t angle_set_num) = 0;

  virtual const double* PsiBndry(uint64_t bndry_map,
                                 int angle_num,
                                 uint64_t cell_local_id,
                                 int face_num,
                                 int fi,
                                 int g,
                                 int gs_ss_begin,
                                 bool surface_source_active) = 0;
  virtual double* ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                             int angle_num,
                                             uint64_t cell_local_id,
                                             int face_num,
                                             int fi,
                                             int gs_ss_begin) = 0;

  virtual ~AngleSet() = default;

protected:
  const size_t num_grps;
  const SPDS& spds_;
  std::shared_ptr<FLUDS> fluds_;
  const std::vector<size_t> angles_;
  std::map<uint64_t, SweepBndryPtr>& ref_boundaries_;
  const size_t ref_group_subset_;

  bool executed = false;
};



} // namespace chi_mesh::sweep_management

#endif // CHI_ANGLESET_H