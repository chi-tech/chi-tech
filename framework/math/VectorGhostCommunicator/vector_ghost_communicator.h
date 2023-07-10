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
protected:
  const uint64_t             local_size_;
  const uint64_t             globl_size_;
  const std::vector<int64_t> ghost_indices_;
  const MPI_Comm             comm_;
  int                        location_id_ = 0;
  int                        process_count_ = 0;
  std::vector<uint64_t>      locI_extents_;
  std::vector<int>           sendcounts_;
  std::vector<int>           senddispls_;
  std::vector<int>           recvcounts_;
  std::vector<int>           recvdispls_;
  std::vector<int64_t>       takers_local_ids_;
  std::map<int64_t,size_t>   ghost_ids_to_recv_map_;

public:
  VectorGhostCommunicator(uint64_t local_size,
                          uint64_t global_size,
                          std::vector<int64_t> ghost_indices,
                          MPI_Comm communicator);

private:
  int FindOwnerPID(uint64_t global_id) const;

public:
  void CommunicateGhostEntries(std::vector<double>& local_vector) const;

  std::vector<double> MakeGhostedVector() const;
  std::vector<double> MakeGhostedVector(
    const std::vector<double>& unghosted_vector) const;
};

}//namespace chi_math

#endif //CHITECH_VECTOR_GHOST_COMMUNICATOR_H
