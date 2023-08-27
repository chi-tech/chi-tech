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
      num_ghosts_(ghost_ids.size()),
      ghost_ids_(ghost_ids),
      ghost_comm_(local_size, global_size, ghost_ids, communicator)
  {}

  GhostedParallelVector(const GhostedParallelVector& other)
    : ParallelVector(other),
      num_ghosts_(other.num_ghosts_),
      ghost_ids_(other.ghost_ids_),
      ghost_comm_(other.ghost_comm_)
  {}


  GhostedParallelVector(GhostedParallelVector&& other)
    : ParallelVector(other),
      num_ghosts_(other.num_ghosts_),
      ghost_ids_(other.ghost_ids_),
      ghost_comm_(other.ghost_comm_)
  {}

  uint64_t NumGhosts() const { return num_ghosts_; }
  uint64_t TotalLocalSize() const override { return local_size_ + num_ghosts_; }

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
  const uint64_t num_ghosts_;
  const std::vector<int64_t> ghost_ids_;
  VectorGhostCommunicator ghost_comm_;


};

}

#endif //GHOSTED_PARALLEL_VECTOR_H
