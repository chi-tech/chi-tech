#ifndef GHOSTED_PARALLEL_VECTOR_H
#define GHOSTED_PARALLEL_VECTOR_H

#include "parallel_vector.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"


namespace chi_math
{


class GhostedParallelVector : public ParallelVector
{
public:
  GhostedParallelVector() = delete;

  /**
   * Initialize the ghosted parallel vector with the given local and global
   * sizes, with the specified global ghost indices.
   */
  GhostedParallelVector(const uint64_t local_size,
                        const uint64_t global_size,
                        const std::vector<int64_t>& ghost_ids,
                        const MPI_Comm communicator = MPI_COMM_WORLD)
    : ParallelVector(local_size, global_size, communicator),
      ghost_comm_(local_size, global_size, ghost_ids, communicator)
  { values_.assign(local_size_ + ghost_comm_.NumGhosts(), 0.0); }

  /// Initialize a ghosted parallel vector with a ghost communicator.
  GhostedParallelVector(const VectorGhostCommunicator& ghost_comm)
    : ParallelVector(ghost_comm.LocalSize(),
                     ghost_comm.GlobalSize(),
                     ghost_comm.Communicator()),
      ghost_comm_(ghost_comm)
  {
    values_.assign(local_size_ + ghost_comm_.NumGhosts(), 0.0);
  }

  GhostedParallelVector(const GhostedParallelVector& other)
    : GhostedParallelVector(other.ghost_comm_)
  {}

  GhostedParallelVector(GhostedParallelVector&& other)
    : GhostedParallelVector(other.ghost_comm_)
  {}

  uint64_t NumGhosts() const { return ghost_comm_.NumGhosts(); }
  uint64_t TotalLocalSize() const override
  { return local_size_ + ghost_comm_.NumGhosts(); }

  const std::vector<int64_t>&
  GhostIndices() const { return ghost_comm_.GhostIndices(); }

  int64_t MapGhostToLocal(const int64_t ghost_id) const
  { return ghost_comm_.MapGhostToLocal(ghost_id); }

  std::vector<double> MakeGhostedLocalVector() const { return values_; }
  
  /**
   * Return the value of the parallel vector for the specified global index.
   *
   * An error is thrown if the global index does not belong to a locally
   * owned, or ghost entry.
   */
  double GetGlobalValue(const int64_t global_id) const;

  /**
   * Communicate the current ghost entries to all other processes to
   * update the locally stored ghost data.
   */
  void CommunicateGhostEntries()
  { ghost_comm_.CommunicateGhostEntries(values_); }

private:
  VectorGhostCommunicator ghost_comm_;


};

}

#endif //GHOSTED_PARALLEL_VECTOR_H
