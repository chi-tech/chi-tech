#ifndef GHOSTED_PARALLEL_VECTOR_H
#define GHOSTED_PARALLEL_VECTOR_H

#include "parallel_vector.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"


namespace chi_math
{


class GhostedParallelVector : public ParallelVector
{
public:
  GhostedParallelVector(const uint64_t local_size,
                        const uint64_t global_size,
                        const std::vector<int64_t>& ghost_indices,
                        const MPI_Comm communicator = MPI_COMM_WORLD);

  uint64_t NumGhosts() const { return ghost_indices_.size(); }

  void CommunicateGhostEntries();

protected:
  const std::vector<int64_t> ghost_indices_;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<int64_t> local_ids_to_send_;
  std::map<int64_t, size_t> ghost_ids_to_recv_map_;

};

}

#endif //GHOSTED_PARALLEL_VECTOR_H
