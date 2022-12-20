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
  const uint64_t             m_local_size;
  const uint64_t             m_globl_size;
  const std::vector<int64_t> m_ghost_indices;
  const MPI_Comm             m_comm;
  int                        m_location_id = 0;
  int                        m_process_count = 0;
  std::vector<uint64_t>      m_locI_extents;
  std::vector<int>           m_sendcounts;
  std::vector<int>           m_senddispls;
  std::vector<int>           m_recvcounts;
  std::vector<int>           m_recvdispls;
  std::vector<int64_t>       m_takers_local_ids;
  std::map<int64_t,size_t>   m_ghost_ids_to_recv_map;

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
