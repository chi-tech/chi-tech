#include "ghosted_parallel_vector.h"

#include "chi_mpi_utils_map_all2all.h"

#include "chi_log.h"
#include "chi_log_exceptions.h"

namespace chi_math
{

GhostedParallelVector::GhostedParallelVector(
    const uint64_t local_size,
    const uint64_t global_size,
    const std::vector<int64_t>& ghost_indices,
    const MPI_Comm communicator)
    : ParallelVector(local_size, global_size, communicator),
      ghost_indices_(ghost_indices)
{
  Reinit(local_size, global_size, ghost_indices);
}


void GhostedParallelVector::Clear()
{
  ParallelVector::Clear();
  ghost_indices_.clear();

  sendcounts_.clear();
  senddispls_.clear();
  recvcounts_.clear();
  recvdispls_.clear();

  local_ids_to_send_.clear();
  ghost_ids_to_recv_map_.clear();
}


void GhostedParallelVector::Reinit(
    const uint64_t local_size,
    const uint64_t global_size,
    const std::vector<int64_t>& ghost_indices)
{
  Clear();
  ParallelVector::Reinit(local_size, global_size);
  ghost_indices_ = ghost_indices;

  // Resize the vector for ghosted entries
  values_.assign(local_size_ + ghost_indices_.size(), 0.0);

  // Construct a mapping between the processes and the ghost indices that
  // need to be sent to the current process by other processes
  std::map<int, std::vector<int64_t>> recv_map;
  for (int64_t ghost_id : ghost_indices)
    recv_map[FindProcessID(ghost_id)].push_back(ghost_id);

  // Now, define an ordering for the ghost indices based off of process.
  // In doing this, all ghost indices that belong to a particular non-local
  // process are ordered contiguously.
  size_t count = 0;
  for (const auto& [pid, ghost_id_list] : recv_map)
    for (int64_t ghost_id : ghost_id_list)
      ghost_ids_to_recv_map_[ghost_id] = count++;

  // Next, the communication structure for receiving data needs to
  // be constructed. This involves determining how much data will be
  // sent by each process, and the starting location of a given process'
  // sent data in the serialized received communication.
  int total_recvcounts = 0;
  recvcounts_.assign(process_count_, 0);
  recvdispls_.assign(process_count_, 0);
  for (const auto& [pid, ghost_id_list] : recv_map)
  {
    recvcounts_[pid] = static_cast<int>(ghost_id_list.size());
    recvdispls_[pid] = total_recvcounts;
    total_recvcounts += ghost_id_list.size();
  }

  // Now, the information about which processes need what information from
  // which other processes is communicated to each other process. This
  // results in each process knowing what it needs to send.
  std::map<int, std::vector<int64_t>> send_map =
      chi_mpi_utils::MapAllToAll(recv_map, MPI_INT64_T, communicator_);

  // With this information, the amount of information that needs
  // to be sent can be determined.
  size_t send_size = 0;
  for (const auto& [pid, ghost_id_list] : send_map)
    send_size += ghost_id_list.size();

  // Next, find the local indices of the data that needs to be sent.
  local_ids_to_send_.reserve(send_size);
  for (const auto& [pid, ghost_ids_list] : send_map)
    for (int64_t ghost_id : ghost_ids_list)
    {
      const auto local_block_beg = extents_[location_id_];
      const auto local_block_end = extents_[location_id_ + 1];

      ChiLogicalErrorIf(
          ghost_id < local_block_beg or ghost_id >= local_block_end,
          std::string(__FUNCTION__) + ": " +
                                    "Problem determining communication pattern. Process " +
                                    std::to_string(pid) + " determined that process " +
                                    std::to_string(location_id_) + " needs to communicate global index " +
                                    std::to_string(ghost_id) + " to it, but this index does not " +
                                    "belong on process " + std::to_string(pid));

      local_ids_to_send_.push_back(ghost_id - local_block_beg);
    }

  // Finally, build the communication structure for data being sent to
  // other processes in a similar manner to that for data being received.
  int total_sendcounts = 0;
  sendcounts_.assign(process_count_, 0);
  senddispls_.assign(process_count_, 0);
  for (const auto& [pid, ghost_id_list] : send_map)
  {
    sendcounts_[pid] = static_cast<int>(ghost_id_list.size());
    senddispls_[pid] = total_sendcounts;
    total_sendcounts += ghost_id_list.size();
  }
}


void GhostedParallelVector::Assemble()
{
  ParallelVector::Assemble();
  CommunicateGhostEntries();
}


void GhostedParallelVector::CommunicateGhostEntries()
{
  // Get the data that needs to be sent
  const size_t send_size = local_ids_to_send_.size();
  std::vector<double> send_data;
  send_data.reserve(send_size);
  for (const int64_t local_id : local_ids_to_send_)
    send_data.push_back(values_[local_id]);

  // Allocate data to be received
  const size_t recv_size = ghost_indices_.size();
  std::vector<double> recv_data(recv_size, 0.0);

  // Communicate the information
  MPI_Alltoallv(send_data.data(),
                sendcounts_.data(),
                senddispls_.data(),
                MPI_DOUBLE,
                recv_data.data(),
                recvcounts_.data(),
                recvdispls_.data(),
                MPI_DOUBLE,
                communicator_);

  // Lastly, populate the local vector with ghost data. All ghost data is
  // appended to the back of the local vector. Using the mapping between
  // ghost indices and the relative ghost index position along with the
  // ordering of the ghost indices, this can be accomplished.
  for (size_t k = 0; k < recv_size; ++k)
    values_[local_size_ + k] =
        recv_data[ghost_ids_to_recv_map_[ghost_indices_[k]]];
}

}

