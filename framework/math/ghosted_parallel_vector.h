#ifndef GHOSTED_PARALLEL_VECTOR_H
#define GHOSTED_PARALLEL_VECTOR_H

#include "parallel_vector.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"


namespace chi_math
{


class GhostedParallelVector : public ParallelVector
{
public:
  /**
   * Initialize a ghosted parallel vector with the given local and global
   * sizes with the given global ghost indices.
   */
  GhostedParallelVector(const uint64_t local_size,
                        const uint64_t global_size,
                        const std::vector<int64_t>& ghost_indices,
                        const MPI_Comm communicator = MPI_COMM_WORLD);

  /// Return the number of ghosts that belong to this process.
  uint64_t NumGhosts() const { return ghost_indices_.size(); }

  /// Clear the ghosted parallel vector, returning it to an unintialized state.
  void Clear() override;

  /**
   * Reinitialize the ghosted parallel vector to the given local and global
   * sizes with the given global ghost indices.
   */
  void Reinit(const uint64_t local_size,
              const uint64_t global_size,
              const std::vector<int64_t>& ghost_indices);

  /**
   * Communicate all operations stored within the operation cache to the
   * corresponding processes that own the respective global indices, then
   * communicate ghost entries.
   *
   * This routine clears the operation cache once completed.
   */
  void Assemble() override;

  /**
   * Communicate the current ghost entries to all other processes to
   * update the locally stored ghost data.
   */
  void CommunicateGhostEntries();

protected:
  std::vector<int64_t> ghost_indices_;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<int64_t> local_ids_to_send_;
  std::map<int64_t, size_t> ghost_ids_to_recv_map_;

};

}

#endif //GHOSTED_PARALLEL_VECTOR_H
