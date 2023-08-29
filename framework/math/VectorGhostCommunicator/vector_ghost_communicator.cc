#include "vector_ghost_communicator.h"

#include "mpi/chi_mpi_utils_map_all2all.h"

#include "chi_log_exceptions.h"

#include <map>
#include <string>
#include <algorithm>


namespace chi_math
{

//######################################################################
VectorGhostCommunicator::VectorGhostCommunicator(
    const uint64_t local_size,
    const uint64_t global_size,
    const std::vector<int64_t>& ghost_ids,
    const MPI_Comm communicator)
    : local_size_(local_size),
      global_size_(global_size),
      ghost_ids_(ghost_ids),
      comm_(communicator)
{
  // Get the process ID and total process count
  MPI_Comm_rank(comm_, &location_id_);
  MPI_Comm_size(comm_, &process_count_);

  // Get the local vector sizes per process
  std::vector<uint64_t> local_sizes(process_count_, 0);
  MPI_Allgather(&local_size_,         //sendbuf
                1, MPI_UINT64_T,      //sendcount + sendtype
                local_sizes.data(),   //recvbuf
                1, MPI_UINT64_T,      //recvcount + recvtype
                comm_);               //communicator

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  extents_.assign(process_count_ + 1, 0);
  for (size_t locJ = 1; locJ < process_count_; ++locJ)
    extents_[locJ] = extents_[locJ - 1] + local_sizes[locJ - 1];
  extents_[process_count_] =
      extents_[process_count_ - 1] + local_sizes.back();

  // Construct a mapping between processes and the ghost indices
  // that belong to them. This information yields what this process
  // needs to receive from other processes.
  std::map<int, std::vector<int64_t>> recv_map;
  for (int64_t ghost_id : ghost_ids)
    recv_map[FindOwnerPID(ghost_id)].push_back(ghost_id);

  // This process will receive data in process-contiguous manner,
  // so a mapping needs to be developed to map each ghost id to
  // its respective ordering in the received data.
  size_t count = 0;
  for (const auto& [pid, gids] : recv_map)
    for (const int64_t gid : gids)
      ghost_to_recv_map_[gid] = count++;

  // Now, the structure of the data being received from communication
  // is developed. This involves determining the amount of data
  // being sent per process and the starting position of the data in
  // the receive buffer per process.
  int total_recvcounts = 0;
  recvcounts_.assign(process_count_, 0);
  recvdispls_.assign(process_count_, 0);
  for (const auto& [pid, gids] : recv_map)
  {
    recvcounts_[pid] = static_cast<int>(gids.size());
    recvdispls_[pid] = total_recvcounts;
    total_recvcounts += gids.size();
  }

  // For communication, each process must also know what it is
  // sending to other processors. If each processor sends each
  // other process the global ids it needs to receive, then each
  // process will know what other processes need from it. The
  // MPI utility MapAllToAll in Chi-Tech accomplishes this task,
  // returning a mapping of processes to the global ids that this
  // process needs to send.
  std::map<int, std::vector<int64_t>> send_map =
      chi_mpi_utils::MapAllToAll(recv_map, MPI_INT64_T, comm_);

  // With this information, the amount of information that needs
  // to be sent can be determined.
  size_t send_size = 0;
  for (const auto& [pid, gids] : send_map)
    send_size += gids.size();

  // Next, the local ids on this process that need to be
  // communicated to other processes can be determined and stored.
  local_ids_to_send_.reserve(send_size);
  for (const auto& [pid, gids] : send_map)
    for (const int64_t gid : gids)
    {
      ChiLogicalErrorIf(
          gid < extents_[location_id_] or
          gid >= extents_[location_id_ + 1],
          std::string(__FUNCTION__) + ": " +
          "Problem determining communication pattern. Process " +
          std::to_string(pid) + " determined that process " +
          std::to_string(location_id_) + " needs to communicate global id " +
          std::to_string(gid) + " to it, but this id is not locally owned.");

      local_ids_to_send_.push_back(gid - extents_[location_id_]);
    }

  // Finally, the communication pattern for the data being sent
  // can be constructed similarly to that for the received data.
  int total_sendcounts = 0;
  sendcounts_.assign(process_count_, 0);
  senddispls_.assign(process_count_, 0);
  for (const auto& [pid, gids] : send_map)
  {
    sendcounts_[pid] = static_cast<int>(gids.size());
    senddispls_[pid] = total_sendcounts;
    total_sendcounts += gids.size();
  }
}


//######################################################################
int64_t VectorGhostCommunicator::
MapGhostToLocal(const int64_t ghost_id) const
{
  ChiInvalidArgumentIf(
      ghost_to_recv_map_.count(ghost_id) == 0,
      "The given ghost id does not belong to this communicator.");

  // Get the position within the ghost id vector of the given ghost id
  const auto k = std::find(ghost_ids_.begin(),
                           ghost_ids_.end(), ghost_id) - ghost_ids_.begin();

  // Local index is local size plus the position in the ghost id vector
  return local_size_ + k;
}


//######################################################################
void VectorGhostCommunicator::
CommunicateGhostEntries(std::vector<double>& ghosted_vector) const
{
  ChiInvalidArgumentIf(
      ghosted_vector.size() != local_size_ + ghost_ids_.size(),
      std::string(__FUNCTION__) + ": Vector size mismatch.");

  // Serialize the data that needs to be sent
  const size_t send_size = local_ids_to_send_.size();
  std::vector<double> send_data;
  send_data.reserve(send_size);
  for (const int64_t local_id : local_ids_to_send_)
    send_data.push_back(ghosted_vector[local_id]);

  // Create serialized storage for the data to be received
  const size_t recv_size = ghost_ids_.size();
  std::vector<double> recv_data(recv_size, 0.0);

  // Communicate the ghost data
  MPI_Alltoallv(send_data.data(),
                sendcounts_.data(),
                senddispls_.data(),
                MPI_DOUBLE,
                recv_data.data(),
                recvcounts_.data(),
                recvdispls_.data(),
                MPI_DOUBLE,
                comm_);

  // Lastly, populate the local vector with ghost data. All ghost data is
  // appended to the back of the local vector. Using the mapping between
  // ghost indices and the relative ghost index position along with the
  // ordering of the ghost indices, this can be accomplished.
  for (size_t k = 0; k < recv_size; ++k)
    ghosted_vector[local_size_ + k] =
        recv_data[ghost_to_recv_map_.at(ghost_ids_[k])];
}


//######################################################################
std::vector<double> VectorGhostCommunicator::MakeGhostedVector() const
{
  const auto ghosted_size = local_size_ + ghost_ids_.size();
  return std::vector<double>(ghosted_size, 0.0);
}


//######################################################################
std::vector<double> VectorGhostCommunicator::
MakeGhostedVector(const std::vector<double>& local_vector) const
{
  ChiInvalidArgumentIf(
      local_vector.size() != local_size_,
      std::string(__FUNCTION__) + ": Incompatible unghosted vector." +
      "unghosted_vector.size() != local_size_");

  // Add ghost indices to the back of the unghosted vector
  std::vector<double> vec = local_vector;
  for (size_t i = 0; i < ghost_ids_.size(); ++i)
    vec.emplace_back(0.0);
  return vec;
}

//###################################################################
int VectorGhostCommunicator::
FindOwnerPID(const int64_t global_id) const
{
  ChiInvalidArgumentIf(
      global_id < 0 or global_id >= global_size_,
      std::string(__FUNCTION__) + ": Invalid global id." +
                                "Global ids must be in [0, global_size_).");

  for (int p = 0; p < process_count_; ++p)
    if (global_id >= extents_[p] and global_id < extents_[p + 1])
      return p;
  return -1;
}

}
