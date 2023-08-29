#ifndef CHITECH_VECTOR_GHOST_COMMUNICATOR_H
#define CHITECH_VECTOR_GHOST_COMMUNICATOR_H

#include <vector>
#include <cstdint>
#include <map>

#include <mpi.h>

namespace chi_math
{

/**Vector with allocation space for ghosts.*/
class VectorGhostCommunicator
{

public:
  VectorGhostCommunicator(const uint64_t local_size,
                          const uint64_t global_size,
                          const std::vector<int64_t>& ghost_ids,
                          const MPI_Comm communicator);

  uint64_t LocalSize() const { return local_size_; }
  uint64_t GlobalSize() const { return global_size_; }

  uint64_t NumGhosts() const { return ghost_ids_.size(); }
  const std::vector<int64_t>& GhostIndices() const { return ghost_ids_; }

  MPI_Comm Communicator() const { return comm_; }

  int64_t MapGhostToLocal(const int64_t ghost_id) const;

  void CommunicateGhostEntries(std::vector<double>& ghosted_vector) const;

  std::vector<double> MakeGhostedVector() const;
  std::vector<double>
  MakeGhostedVector(const std::vector<double>& local_vector) const;

private:
  int FindOwnerPID(const int64_t global_id) const;

protected:
  const uint64_t local_size_;
  const uint64_t global_size_;
  const std::vector<int64_t> ghost_ids_;
  const MPI_Comm comm_;

  int location_id_;
  int process_count_;
  std::vector<uint64_t> extents_;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<int64_t> local_ids_to_send_;
  std::map<int64_t, size_t> ghost_to_recv_map_;

};

}//namespace chi_math

#endif //CHITECH_VECTOR_GHOST_COMMUNICATOR_H
